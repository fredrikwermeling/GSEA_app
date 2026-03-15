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

        // Gene set browser state
        this.selectedGeneSets = new Set();     // Set of gene set names for custom selection
        this.useCustomSelection = false;        // false = collection checkboxes; true = per-set selection
        this._gsbAllItems = [];                 // master sorted list built when modal opens
        this._gsbFlatList = [];                 // filtered flat list (headers + items)
        this._gsbVisualRows = [];               // visual rows for virtual scroll (multi-column)
        this._gsbColumns = 3;                   // number of columns for gene set items
        this._gsbRowHeight = 32;
        this._gsbBufferRows = 12;

        // Settings
        this.settings = {
            permutations: 1000,
            minSize: 15,
            maxSize: 500,
            weightP: 1,
            topN: 20,
            fdrDisplayThreshold: 0.25,
            colorScale: 'YlOrRd',
            exportScale: 4,
            transparentBg: false,
            // Figure customization
            fontFamily: 'Open Sans',
            fontSize: 12,
            esLineColor: '#15a04a',
            esLineWidth: 2.5,
            showStatsBox: true,
            showZeroCross: false,
            showESIndicator: true,
            showCorrelationLabels: true,
            showPanelBorders: true,
            showHitMarkers: true,
            showMetricFill: true,
            positiveColor: '#dc2626',
            negativeColor: '#2563eb',
            overlapColorScheme: 'green',
            dataType: 'expression',      // 'expression' | 'crispr'
            // Height controls
            esPlotWidth: 500,
            esPlotHeight: 500,
            bubblePlotHeight: 0,         // 0 = auto
            rankedPlotHeight: 250,
            overlapPlotHeight: 0,        // 0 = auto
            // ES panel proportions (%)
            esPanelES: 68,
            esPanelHits: 10,
            esPanelMetric: 20
        };

        // Per-text-element font settings
        this.textElements = {
            bubble: [
                { key: 'title', label: 'Title', editable: true, defaultText: 'Enrichment Overview', defaultSize: 14 },
                { key: 'xAxisLabel', label: 'X-axis label', editable: true, defaultText: 'Normalized Enrichment Score (NES)', defaultSize: 12 },
                { key: 'yTickFont', label: 'Y-axis labels', editable: false, defaultSize: 11 },
                { key: 'colorbarTitle', label: 'Colorbar title', editable: true, defaultText: 'FDR', defaultSize: 11 },
                { key: 'colorbarTickFont', label: 'Colorbar ticks', editable: false, defaultSize: 10 },
                { key: 'sizeAnnotation', label: 'Size annotation', editable: false, defaultSize: 10 }
            ],
            ranked: [
                { key: 'title', label: 'Title', editable: true, defaultText: 'Ranked List Metric', defaultSize: 14 },
                { key: 'xAxisLabel', label: 'X-axis label', editable: true, defaultText: 'Rank in Ordered Dataset', defaultSize: 11 },
                { key: 'yAxisLabel', label: 'Y-axis label', editable: true, defaultText: 'Ranked list metric', defaultSize: 11 },
                { key: 'xTickFont', label: 'X-axis ticks', editable: false, defaultSize: 10 },
                { key: 'yTickFont', label: 'Y-axis ticks', editable: false, defaultSize: 10 },
                { key: 'posLabel', label: 'Positive label', editable: true, defaultText: 'Positively correlated', defaultSize: 9 },
                { key: 'negLabel', label: 'Negative label', editable: true, defaultText: 'Negatively correlated', defaultSize: 9 }
            ],
            es: [
                { key: 'title', label: 'Title', editable: true, defaultSize: 16 },
                { key: 'esYLabel', label: 'ES Y-axis label', editable: true, defaultText: 'Enrichment score (ES)', defaultSize: 14 },
                { key: 'hitMarkersLabel', label: 'Hit markers label', editable: true, defaultText: 'Gene hits', defaultSize: 10 },
                { key: 'metricYLabel', label: 'Metric Y-axis label', editable: true, defaultText: '', defaultSize: 13 },
                { key: 'xAxisLabel', label: 'X-axis label', editable: true, defaultText: 'Rank in Ordered Dataset', defaultSize: 14 },
                { key: 'statsBox', label: 'Stats box', editable: false, defaultSize: 13, defaultFamily: 'Roboto Mono' },
                { key: 'posLabel', label: 'Positive label', editable: true, defaultText: 'Positively correlated', defaultSize: 12 },
                { key: 'negLabel', label: 'Negative label', editable: true, defaultText: 'Negatively correlated', defaultSize: 12 },
                { key: 'esTickFont', label: 'ES axis ticks', editable: false, defaultSize: 11 },
                { key: 'tickFont', label: 'Metric axis ticks', editable: false, defaultSize: 11 }
            ],
            overlap: [
                { key: 'title', label: 'Title', editable: true, defaultText: 'Gene Set Overlap', defaultSize: 14 },
                { key: 'colorbarTitle', label: 'Colorbar title', editable: true, defaultText: 'Jaccard Index', defaultSize: 11 },
                { key: 'tickFont', label: 'Axis labels', editable: false, defaultSize: 10 }
            ]
        };
        this.textFonts = {};
        this._initTextFonts();

        // Stats box corner cache (preserves position across re-renders)
        this._statsBoxCornerCache = {};
        this._textScaleOffset = { bubble: 0, ranked: 0, es: 0, overlap: 0 };

        // Overlap-based redundancy filtering
        this._overlapCache = null;       // { pairwise: Map, setGenes: {} }
        this._hiddenSets = new Set();    // gene set names hidden by user or overlap filter
        this._pinnedBubbleSets = new Set(); // gene sets pinned to always show in bubble plot
        this._overlapFilterThreshold = 0; // 0 = off, 50/70/80 = % Jaccard overlap

        // Table sort state
        this.sortCol = 'fdr';
        this.sortAsc = true;

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

        // R results upload
        document.getElementById('rResultsInput').addEventListener('change', (e) => {
            if (e.target.files[0]) this.loadRResults(e.target.files[0]);
        });

        // Column selection
        document.getElementById('geneColumn').addEventListener('change', () => this.checkReady());
        document.getElementById('metricColumn').addEventListener('change', () => this.checkReady());

        // Data type selector (sidebar + inline sync)
        document.getElementById('dataType').addEventListener('change', (e) => {
            this.settings.dataType = e.target.value;
            document.getElementById('dataTypeInline').value = e.target.value;
            const sel = document.getElementById('geneSetSelector');
            if (sel && sel.value) this.renderGeneSetInfo(sel.value);
            if (this.results) this.updateSettings();
        });
        document.getElementById('dataTypeInline').addEventListener('change', (e) => {
            this.settings.dataType = e.target.value;
            document.getElementById('dataType').value = e.target.value;
            const sel = document.getElementById('geneSetSelector');
            if (sel && sel.value) this.renderGeneSetInfo(sel.value);
            if (this.results) this.updateSettings();
        });

        // Gene set collection checkboxes
        ['checkHallmark', 'checkC2kegg', 'checkC2reactome', 'checkC2wp', 'checkC2biocarta', 'checkC2pid', 'checkC2cgp', 'checkC3', 'checkC5bp', 'checkC5cc', 'checkC5mf', 'checkC5hpo', 'checkC6', 'checkC7', 'checkC8'].forEach(id => {
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

        // Gene Set Browser
        document.getElementById('openGeneSetBrowser').addEventListener('click', (e) => {
            e.preventDefault();
            this.openGeneSetBrowser();
        });
        document.getElementById('gsbClose').addEventListener('click', () => this.closeGeneSetBrowser());
        document.getElementById('gsbBackdrop').addEventListener('click', () => this.closeGeneSetBrowser());
        document.getElementById('gsbCancelBtn').addEventListener('click', () => this.closeGeneSetBrowser());
        document.getElementById('gsbApplyBtn').addEventListener('click', () => this.applyGeneSetSelection());
        document.getElementById('gsbCollectionFilter').addEventListener('change', () => this.filterGeneSetBrowser());
        document.getElementById('gsbTierFilter').addEventListener('change', () => {
            this._updateCollectionFilterForTier();
            this.filterGeneSetBrowser();
        });
        document.getElementById('gsbJaccardFilter').addEventListener('change', () => this.filterGeneSetBrowser());
        document.getElementById('gsbSelectAll').addEventListener('click', () => this.gsbSelectAllVisible());
        document.getElementById('gsbDeselectAll').addEventListener('click', () => this.gsbDeselectAll());
        document.getElementById('gsbDownloadBtn').addEventListener('click', () => this.gsbDownloadSelectedCSV());

        // Gene Set Browser: debounced search
        let gsbSearchTimer = null;
        document.getElementById('gsbSearch').addEventListener('input', () => {
            clearTimeout(gsbSearchTimer);
            gsbSearchTimer = setTimeout(() => this.filterGeneSetBrowser(), 150);
        });

        // Gene Set Browser: virtual scroll with requestAnimationFrame
        let gsbRafPending = false;
        document.getElementById('gsbListContainer').addEventListener('scroll', () => {
            if (!gsbRafPending) {
                gsbRafPending = true;
                requestAnimationFrame(() => {
                    this.gsbRenderVisibleRows();
                    gsbRafPending = false;
                });
            }
        });

        // Gene Set Browser: event delegation for row interactions
        const gsbRowsEl = document.getElementById('gsbRows');
        gsbRowsEl.addEventListener('click', (e) => {
            const row = e.target.closest('.gsb-row');
            if (!row) return;
            const checkbox = row.querySelector('input[type="checkbox"]');
            if (!checkbox) return;
            const name = checkbox.dataset.name;
            if (e.target !== checkbox) {
                checkbox.checked = !checkbox.checked;
            }
            if (checkbox.checked) {
                this.selectedGeneSets.add(name);
                row.classList.add('selected');
            } else {
                this.selectedGeneSets.delete(name);
                row.classList.remove('selected');
            }
            this.updateGsbSelectionCount();
        });
        gsbRowsEl.addEventListener('mouseover', (e) => {
            const row = e.target.closest('.gsb-row');
            if (!row) return;
            const name = row.dataset.name;
            if (!name || !this._gsbItemMap) return;
            const item = this._gsbItemMap.get(name);
            if (item) this.gsbShowDetail(item);
        });

        // Close gene set browser on Escape
        document.addEventListener('keydown', (e) => {
            if (e.key === 'Escape' && document.getElementById('gsbModal').classList.contains('open')) {
                this.closeGeneSetBrowser();
            }
        });

        // Run / Cancel / Smart Run
        document.getElementById('runBtn').addEventListener('click', () => this.runGSEA());
        document.getElementById('smartRunBtn').addEventListener('click', () => this.runSmartGSEA());
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

        // Search filter for gene set selector
        this._initGeneSetSearch();

        // Enrichment Plot tab filters — re-populate dropdown when changed
        ['esFdrFilter', 'esPvalFilter', 'esDirectionFilter', 'esRedundancyFilter'].forEach(id => {
            document.getElementById(id).addEventListener('change', () => {
                // Clear search filter when changing FDR/pval/direction
                const searchInput = document.getElementById('geneSetSearchInput');
                if (searchInput) searchInput.value = '';
                this.populateGeneSetSelector();
                this._initGeneSetSearch();
            });
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
        // Live-update overlap heatmap on parameter change (debounced)
        let overlapDebounce = null;
        const overlapLiveUpdate = () => {
            clearTimeout(overlapDebounce);
            overlapDebounce = setTimeout(() => { if (this.results) this.renderOverlapHeatmap(); }, 300);
        };
        document.getElementById('overlapMaxSets').addEventListener('input', overlapLiveUpdate);
        ['overlapFdrFilter', 'overlapPvalFilter', 'overlapDirectionFilter', 'overlapRedundancyFilter'].forEach(id => {
            const el = document.getElementById(id);
            if (el) el.addEventListener('change', overlapLiveUpdate);
        });

        // Overlap → Results Table: show only the sets visible in the overlap heatmap
        document.getElementById('overlapToTableBtn').addEventListener('click', () => {
            if (!this._overlapVisibleSets) return;
            // Hide all sets, then show only those in the overlap heatmap
            this._hiddenSets.clear();
            for (const r of this.results) {
                if (!this._overlapVisibleSets.has(r.name)) {
                    this._hiddenSets.add(r.name);
                }
            }
            this.renderBubblePlot();
            this.filterAndRenderTable();
            // Switch to results table tab
            document.querySelector('[data-tab="table"]').click();
        });

        // Table filter
        document.getElementById('tableFilter').addEventListener('input', () => this.filterAndRenderTable());
        document.getElementById('fdrFilter').addEventListener('change', () => this.filterAndRenderTable());
        document.getElementById('pvalueFilter').addEventListener('change', () => this.filterAndRenderTable());
        document.getElementById('directionFilter').addEventListener('change', () => this.filterAndRenderTable());
        document.getElementById('tableOverlapFilter').addEventListener('change', () => this.filterAndRenderTable());

        // Gene set filter popup backdrop
        document.getElementById('geneSetFilterBackdrop').addEventListener('click', () => {
            document.getElementById('geneSetFilterPopup').style.display = 'none';
            document.getElementById('geneSetFilterBackdrop').style.display = 'none';
        });

        // Table column sort
        document.querySelectorAll('#resultsTable thead th').forEach(th => {
            th.addEventListener('click', () => this.sortTable(th.dataset.col));
        });

        // Download CSV
        document.getElementById('downloadCSVBtn').addEventListener('click', () => this.downloadCSV());

        // Re-run filtered gene sets (results table button)
        document.getElementById('rerunFilteredBtn').addEventListener('click', () => this.rerunFiltered());

        // Re-run significant hits (Card 4 button)
        document.getElementById('rerunSignificantBtn').addEventListener('click', () => this.rerunSignificant());
        document.getElementById('expandHitsBtn').addEventListener('click', () => this.expandAroundHits());
        document.getElementById('rerunFdrThreshold').addEventListener('change', () => this.updateRerunHitCount());

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

        // Interpret Guide modal
        const interpretLink = document.getElementById('interpretGuideLink');
        const interpretPopup = document.getElementById('interpretGuidePopup');
        const interpretBackdrop = document.getElementById('interpretGuideBackdrop');
        const interpretClose = document.getElementById('interpretGuideClose');
        interpretLink.addEventListener('click', () => {
            interpretPopup.classList.add('open');
            interpretBackdrop.classList.add('open');
        });
        const closeInterpret = () => {
            interpretPopup.classList.remove('open');
            interpretBackdrop.classList.remove('open');
        };
        interpretClose.addEventListener('click', closeInterpret);
        interpretBackdrop.addEventListener('click', closeInterpret);

        // Changelog modal — close on backdrop click
        const changelogModal = document.getElementById('changelogModal');
        if (changelogModal) {
            changelogModal.addEventListener('click', (e) => {
                if (e.target === changelogModal) changelogModal.style.display = 'none';
            });
        }

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

        // Overview toolbar controls — re-render lollipop on change
        ['overviewFdrFilter', 'overviewPvalFilter', 'overviewDirectionFilter', 'overviewTopN', 'overviewOverlapFilter'].forEach(id => {
            const el = document.getElementById(id);
            if (el) el.addEventListener('change', () => {
                // Sync topN and FDR with gear popup
                if (id === 'overviewTopN') document.getElementById('topN').value = el.value;
                if (id === 'overviewFdrFilter') document.getElementById('fdrDisplayThreshold').value = el.value;
                if (id === 'overviewPvalFilter') document.getElementById('pvalDisplayThreshold').value = el.value;
                if (id === 'overviewOverlapFilter') this._overlapFilterThreshold = parseInt(el.value);
                this.renderBubblePlot();
            });
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

        // Invert ranking & re-run
        document.getElementById('invertBtn').addEventListener('click', () => this.invertAndRerun());

        // Floating settings panel open/close
        document.getElementById('openSettingsBtn').addEventListener('click', () => {
            const panel = document.getElementById('settingsPanel');
            panel.classList.toggle('open');
        });
        document.getElementById('closeSettingsBtn').addEventListener('click', () => {
            document.getElementById('settingsPanel').classList.remove('open');
        });

        // Close popups when clicking outside
        document.addEventListener('click', (e) => {
            // Graph settings popups
            if (!e.target.closest('.card-settings-btn') && !e.target.closest('.graph-settings-popup')) {
                document.querySelectorAll('.graph-settings-popup.open').forEach(p => p.classList.remove('open'));
            }
            // Text settings panel
            const tsp = document.getElementById('textSettingsPanel');
            if (tsp && tsp.style.display !== 'none' && !e.target.closest('.text-settings-panel') && !e.target.closest('button[onclick*="openTextSettings"]')) {
                this.closeTextSettings();
            }
            // Global settings panel
            const sp = document.getElementById('settingsPanel');
            if (sp && sp.classList.contains('open') && !e.target.closest('#settingsPanel') && !e.target.closest('#openSettingsBtn')) {
                sp.classList.remove('open');
            }
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
        // Reset custom selection when user toggles collection checkboxes
        this.useCustomSelection = false;
        this.selectedGeneSets.clear();
        document.getElementById('checkCustomSelection').checked = false;
        document.getElementById('customSelectionCount').textContent = '';

        const collections = {
            checkHallmark: { id: 'hallmark', file: 'h.all.v2023.2.Hs.json' },
            checkC2kegg: { id: 'c2kegg', file: 'c2.cp.kegg.v2023.2.Hs.json' },
            checkC2reactome: { id: 'c2reactome', file: 'c2.cp.reactome.v2023.2.Hs.json' },
            checkC2wp: { id: 'c2wp', file: 'c2.cp.wikipathways.v2023.2.Hs.json' },
            checkC2biocarta: { id: 'c2biocarta', file: 'c2.cp.biocarta.v2023.2.Hs.json' },
            checkC2pid: { id: 'c2pid', file: 'c2.cp.pid.v2023.2.Hs.json' },
            checkC2cgp: { id: 'c2cgp', file: 'c2.cgp.v2023.2.Hs.json' },
            checkC3: { id: 'c3', file: 'c3.all.v2023.2.Hs.json' },
            checkC5bp: { id: 'c5bp', file: 'c5.go.bp.v2023.2.Hs.json' },
            checkC5cc: { id: 'c5cc', file: 'c5.go.cc.v2023.2.Hs.json' },
            checkC5mf: { id: 'c5mf', file: 'c5.go.mf.v2023.2.Hs.json' },
            checkC5hpo: { id: 'c5hpo', file: 'c5.hpo.v2023.2.Hs.json' },
            checkC6: { id: 'c6', file: 'c6.all.v2023.2.Hs.json' },
            checkC7: { id: 'c7', file: 'c7.all.v2023.2.Hs.json' },
            checkC8: { id: 'c8', file: 'c8.all.v2023.2.Hs.json' }
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
        // Custom selection mode
        if (this.useCustomSelection) {
            const n = this.selectedGeneSets.size;
            this.showStatus('geneSetStatus', 'info',
                `Custom selection: ${n.toLocaleString()} gene sets selected — click "Browse" to modify`);
            return;
        }

        const parts = [];
        let total = 0;
        const collections = {
            checkHallmark: { id: 'hallmark', label: 'Hallmark' },
            checkC2kegg: { id: 'c2kegg', label: 'KEGG' },
            checkC2reactome: { id: 'c2reactome', label: 'Reactome' },
            checkC2wp: { id: 'c2wp', label: 'WikiPathways' },
            checkC2biocarta: { id: 'c2biocarta', label: 'BioCarta' },
            checkC2pid: { id: 'c2pid', label: 'PID' },
            checkC2cgp: { id: 'c2cgp', label: 'CGP' },
            checkC3: { id: 'c3', label: 'C3' },
            checkC5bp: { id: 'c5bp', label: 'GO:BP' },
            checkC5cc: { id: 'c5cc', label: 'GO:CC' },
            checkC5mf: { id: 'c5mf', label: 'GO:MF' },
            checkC5hpo: { id: 'c5hpo', label: 'HPO' },
            checkC6: { id: 'c6', label: 'C6' },
            checkC7: { id: 'c7', label: 'C7' },
            checkC8: { id: 'c8', label: 'C8' }
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
            if (total > 2000) {
                this.showStatus('geneSetStatus', 'warning',
                    `${parts.join(' | ')} — ${total.toLocaleString()} total gene sets. ` +
                    `Tip: start with Hallmark (50 sets) for a quick overview. For large runs, you can lower permutations for faster screening, but note that borderline hits (FDR ~0.2–0.3) may be missed — re-run interesting results with ≥1000 permutations.`);
            } else {
                this.showStatus('geneSetStatus', 'info',
                    `${parts.join(' | ')} — ${total.toLocaleString()} total gene sets`);
            }
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
        const smartBtn = document.getElementById('smartRunBtn');
        smartBtn.disabled = !(hasData && hasGeneSets);
        // Only show Smart Run when >1000 gene sets selected
        const geneSets = hasGeneSets ? this.getActiveGeneSets() : null;
        const nSets = geneSets ? Object.keys(geneSets).length : 0;
        smartBtn.style.display = (hasData && hasGeneSets && nSets > 1000) ? '' : 'none';
    }

    getActiveGeneSets() {
        // Custom per-set selection mode
        if (this.useCustomSelection && this.selectedGeneSets.size > 0) {
            const combined = {};
            const allData = { hallmark: this.geneSets['hallmark'], c2kegg: this.geneSets['c2kegg'], c2reactome: this.geneSets['c2reactome'], c2wp: this.geneSets['c2wp'], c2biocarta: this.geneSets['c2biocarta'], c2pid: this.geneSets['c2pid'], c2cgp: this.geneSets['c2cgp'], c3: this.geneSets['c3'], c5bp: this.geneSets['c5bp'], c5cc: this.geneSets['c5cc'], c5mf: this.geneSets['c5mf'], c5hpo: this.geneSets['c5hpo'], c6: this.geneSets['c6'], c7: this.geneSets['c7'], c8: this.geneSets['c8'], custom: this.customGeneSets };
            for (const collData of Object.values(allData)) {
                if (!collData) continue;
                for (const [name, genes] of Object.entries(collData)) {
                    if (this.selectedGeneSets.has(name)) {
                        combined[name] = genes;
                    }
                }
            }
            return Object.keys(combined).length > 0 ? combined : null;
        }

        // Default: whole-collection mode (original behavior)
        const all = {};
        const collections = {
            checkHallmark: 'hallmark',
            checkC2kegg: 'c2kegg',
            checkC2reactome: 'c2reactome',
            checkC2wp: 'c2wp',
            checkC2biocarta: 'c2biocarta',
            checkC2pid: 'c2pid',
            checkC2cgp: 'c2cgp',
            checkC3: 'c3',
            checkC5bp: 'c5bp',
            checkC5cc: 'c5cc',
            checkC5mf: 'c5mf',
            checkC5hpo: 'c5hpo',
            checkC6: 'c6',
            checkC7: 'c7',
            checkC8: 'c8'
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
    // Invert Ranking — negate all metrics and re-run GSEA
    // --------------------------------------------------------
    invertAndRerun() {
        if (!this.rankedList) return;
        // Negate all metrics (this reverses the ranking order)
        this.rankedList.metrics = this.rankedList.metrics.map(m => -m);
        // Reverse the arrays so they stay sorted descending
        this.rankedList.genes.reverse();
        this.rankedList.metrics.reverse();

        // Re-run GSEA with inverted ranking
        const geneSets = this.getActiveGeneSets();
        if (!geneSets) return;
        this.readSettings();

        document.getElementById('runBtn').style.display = 'none';
        document.getElementById('rerunSection').style.display = 'none';
        document.getElementById('cancelBtn').style.display = '';
        const container = document.getElementById('progressContainer');
        container.classList.add('active');
        document.getElementById('progressBar').style.width = '0%';
        document.getElementById('progressText').textContent = 'Re-running with inverted ranking...';
        this.hideStatus('runStatus');

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

    // --------------------------------------------------------
    // GSEA Execution
    // --------------------------------------------------------
    async runGSEA() {
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
        if (nSets > 2000) {
            const action = await this._showRunWarningDialog(nSets, geneSets);
            if (action === 'cancel' || action === 'rscript') {
                container.classList.remove('active');
                document.getElementById('runBtn').style.display = '';
                document.getElementById('cancelBtn').style.display = 'none';
                return;
            }
            // action === 'run' — continue with (possibly updated) permutations
            document.getElementById('progressText').textContent =
                `Processing ${nSets.toLocaleString()} gene sets — this may take several minutes...`;
        }

        this.analysisDate = new Date();
        this._allRunGeneSets = geneSets; // Store for "Expand around hits"

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
        this._smartRunCancelled = true;
        this.createWorker();
        document.getElementById('runBtn').style.display = '';
        document.getElementById('cancelBtn').style.display = 'none';
        document.getElementById('progressContainer').classList.remove('active');
        this.showStatus('runStatus', 'warning', 'Analysis cancelled.');
    }

    // --------------------------------------------------------
    // Run Warning Dialog (replaces native confirm)
    // --------------------------------------------------------
    _showRunWarningDialog(nSets, geneSets) {
        return new Promise((resolve) => {
            const dialog = document.getElementById('runWarningDialog');
            const backdrop = document.getElementById('runWarningBackdrop');
            const permInput = document.getElementById('runWarningPermInput');

            permInput.value = this.settings.permutations;

            document.getElementById('runWarningInfo').textContent =
                `You are about to analyze ${nSets.toLocaleString()} gene sets with ${this.settings.permutations.toLocaleString()} permutations.`;

            const permTip = this.settings.permutations <= 200
                ? `You're using ${this.settings.permutations} permutations (good for screening). After results, use "Re-run filtered" with \u22651000 for publication-quality statistics.`
                : `Lower permutations (e.g. 100\u2013200) for faster screening, then use "Re-run filtered" with \u22651000 permutations.`;

            document.getElementById('runWarningTips').innerHTML =
                `<li>Start with Hallmark (50 sets) for a quick, interpretable overview.</li>
                 <li>${permTip}</li>
                 <li>Consider "Smart Run" for iterative analysis of large collections.</li>
                 <li>Or download an <b>R script</b> to run fgsea locally \u2014 handles any size without browser limits.</li>`;

            const updateEstimate = () => {
                const p = parseInt(permInput.value) || 1000;
                const est = Math.max(1, Math.round(nSets * p / 500000));
                document.getElementById('runWarningEstimate').textContent =
                    `Estimated time: ~${est} minute${est > 1 ? 's' : ''}`;
            };
            updateEstimate();
            permInput.oninput = updateEstimate;

            dialog.style.display = 'block';
            backdrop.style.display = 'block';

            const close = (action) => {
                dialog.style.display = 'none';
                backdrop.style.display = 'none';
                permInput.oninput = null;
                resolve(action);
            };

            const wireBtn = (id, action) => {
                const btn = document.getElementById(id);
                const newBtn = btn.cloneNode(true);
                btn.parentNode.replaceChild(newBtn, btn);
                newBtn.addEventListener('click', () => {
                    if (action === 'run') {
                        this.settings.permutations = parseInt(permInput.value) || 1000;
                        document.getElementById('permutations').value = this.settings.permutations;
                    }
                    if (action === 'rscript') {
                        this.downloadFgseaScript(geneSets);
                    }
                    close(action);
                });
            };
            wireBtn('runWarningCancelBtn', 'cancel');
            wireBtn('runWarningClose', 'cancel');
            wireBtn('runWarningRunBtn', 'run');
            wireBtn('runWarningRScriptBtn', 'rscript');

            backdrop.onclick = () => close('cancel');
        });
    }

    _showSmartRunDialog(nTotal) {
        return new Promise((resolve) => {
            const dialog = document.getElementById('runWarningDialog');
            const backdrop = document.getElementById('runWarningBackdrop');
            const permInput = document.getElementById('runWarningPermInput');

            // Smart Run uses fixed permutations: 100 for screening, 1000 for final
            permInput.value = 1000;
            permInput.parentElement.style.display = 'none';

            document.getElementById('runWarningInfo').innerHTML =
                `<b>Smart Run</b> will analyze <b>${nTotal.toLocaleString()}</b> gene sets in 3 phases:` +
                `<ol style="margin: 6px 0 0 16px; padding: 0; font-size: 0.9em; line-height: 1.5;">` +
                `<li><b>Screen</b> — select diverse (non-redundant) subset, test with 100 permutations</li>` +
                `<li><b>Expand</b> — test sets related to any hits from Phase 1</li>` +
                `<li><b>Refine</b> — re-run significant hits with 1,000 permutations for precise FDR</li></ol>` +
                `<div style="margin-top: 6px; font-size: 0.85em; color: var(--gray-500);">Redundant gene sets (Jaccard overlap > 10%) are skipped since they would produce near-identical results. This is 2–7× faster than testing all sets. Use "Run GSEA" if you want to test every set.</div>`;

            document.getElementById('runWarningTips').innerHTML =
                `<li>You can cancel between phases without losing earlier results.</li>`;

            const est = Math.max(1, Math.round(nTotal * 0.3 * 1000 / 500000));
            document.getElementById('runWarningEstimate').textContent =
                `Estimated time: ~${est} minute${est > 1 ? 's' : ''}`;

            // Hide R script button for smart run
            document.getElementById('runWarningRScriptBtn').style.display = 'none';
            document.getElementById('runWarningRunBtn').textContent = 'Start Smart Run';

            dialog.style.display = 'block';
            backdrop.style.display = 'block';

            const close = (action) => {
                dialog.style.display = 'none';
                backdrop.style.display = 'none';
                permInput.oninput = null;
                // Restore buttons and permutation input visibility
                permInput.parentElement.style.display = '';
                document.getElementById('runWarningRScriptBtn').style.display = '';
                document.getElementById('runWarningRunBtn').textContent = 'Run in Browser';
                resolve(action);
            };

            const wireBtn = (id, action) => {
                const btn = document.getElementById(id);
                const newBtn = btn.cloneNode(true);
                btn.parentNode.replaceChild(newBtn, btn);
                newBtn.addEventListener('click', () => {
                    if (action === 'run') {
                        this.settings.permutations = parseInt(permInput.value) || 1000;
                        document.getElementById('permutations').value = this.settings.permutations;
                    }
                    close(action);
                });
            };
            wireBtn('runWarningCancelBtn', 'cancel');
            wireBtn('runWarningClose', 'cancel');
            wireBtn('runWarningRunBtn', 'run');

            backdrop.onclick = () => close('cancel');
        });
    }

    // --------------------------------------------------------
    // fgsea R Script Generation
    // --------------------------------------------------------
    downloadFgseaScript(geneSets) {
        this.readSettings();
        const genes = this.rankedList.genes;
        const metrics = this.rankedList.metrics;
        const { permutations, minSize, maxSize } = this.settings;

        // Determine which standard MSigDB collections are selected
        const msigdbMap = {
            checkHallmark: { category: 'H', subcategory: null, collId: 'hallmark', label: 'Hallmark' },
            checkC2kegg: { category: 'C2', subcategory: 'CP:KEGG', collId: 'c2kegg', label: 'C2 KEGG' },
            checkC2reactome: { category: 'C2', subcategory: 'CP:REACTOME', collId: 'c2reactome', label: 'C2 Reactome' },
            checkC2wp: { category: 'C2', subcategory: 'CP:WIKIPATHWAYS', collId: 'c2wp', label: 'C2 WikiPathways' },
            checkC2biocarta: { category: 'C2', subcategory: 'CP:BIOCARTA', collId: 'c2biocarta', label: 'C2 BioCarta' },
            checkC2pid: { category: 'C2', subcategory: 'CP:PID', collId: 'c2pid', label: 'C2 PID' },
            checkC2cgp: { category: 'C2', subcategory: 'CGP', collId: 'c2cgp', label: 'C2 CGP' },
            checkC3: { category: 'C3', subcategory: null, collId: 'c3', label: 'C3 (Regulatory)' },
            checkC5bp: { category: 'C5', subcategory: 'GO:BP', collId: 'c5bp', label: 'C5 GO:BP' },
            checkC5cc: { category: 'C5', subcategory: 'GO:CC', collId: 'c5cc', label: 'C5 GO:CC' },
            checkC5mf: { category: 'C5', subcategory: 'GO:MF', collId: 'c5mf', label: 'C5 GO:MF' },
            checkC5hpo: { category: 'C5', subcategory: 'HPO', collId: 'c5hpo', label: 'C5 HPO' },
            checkC6: { category: 'C6', subcategory: null, collId: 'c6', label: 'C6 (Oncogenic)' },
            checkC7: { category: 'C7', subcategory: null, collId: 'c7', label: 'C7 (Immunologic)' },
            checkC8: { category: 'C8', subcategory: null, collId: 'c8', label: 'C8 (Cell Type)' }
        };

        const msigdbCollections = [];
        const msigdbSetNames = new Set();
        for (const [checkId, info] of Object.entries(msigdbMap)) {
            const cb = document.getElementById(checkId);
            if (cb && cb.checked && this.geneSets[info.collId]) {
                msigdbCollections.push(info);
                // Track which gene set names come from standard collections
                if (this.geneSets[info.collId]) {
                    for (const name of Object.keys(this.geneSets[info.collId])) {
                        msigdbSetNames.add(name);
                    }
                }
            }
        }

        // Check for custom gene sets (not from msigdbr)
        const customSets = {};
        for (const [name, setGenes] of Object.entries(geneSets)) {
            if (!msigdbSetNames.has(name)) {
                customSets[name] = setGenes;
            }
        }
        const hasCustomSets = Object.keys(customSets).length > 0;
        const useMsigdbr = msigdbCollections.length > 0;

        // Build ranked list
        let script = `# ============================================================
# fgsea Analysis Script — Generated by Enrich
# ${new Date().toISOString().split('T')[0]}
# https://fredrikwermeling.github.io/GSEA/
# ============================================================
#
# This script runs GSEA using the fgsea R package on your data.
# Results are saved as JSON for re-import into Enrich.
#${useMsigdbr ? `
# Gene sets are fetched from MSigDB via the msigdbr package,
# keeping this script small and reproducible.` : `
# Gene sets are embedded inline (dictionary-encoded).`}
#
# Install dependencies (run once):
#   install.packages("BiocManager")
#   BiocManager::install("fgsea")
#   install.packages("jsonlite")${useMsigdbr ? '\n#   install.packages("msigdbr")' : ''}
# ============================================================

library(fgsea)
library(jsonlite)${useMsigdbr ? '\nlibrary(msigdbr)' : ''}

# --- Ranked gene list (${genes.length} genes) ---
ranked_stats <- c(
`;
        const entries = [];
        for (let i = 0; i < genes.length; i++) {
            entries.push(`    "${this._rEscape(genes[i].toUpperCase())}" = ${metrics[i]}`);
        }
        script += entries.join(',\n') + '\n)\n\n';

        // Gene sets section
        if (useMsigdbr) {
            script += `# --- Gene sets from MSigDB (via msigdbr) ---\n`;
            script += `cat("Fetching gene sets from MSigDB...\\n")\n`;
            script += `gene_sets <- list()\n`;
            for (const coll of msigdbCollections) {
                script += `{\n`;
                script += `    cat("  Loading ${coll.label}...\\n")\n`;
                script += `    m <- msigdbr(species = "Homo sapiens", collection = "${coll.category}"${coll.subcategory ? `, subcollection = "${coll.subcategory}"` : ''})\n`;
                script += `    sets <- split(m$gene_symbol, m$gs_name)\n`;
                script += `    gene_sets <- c(gene_sets, sets)\n`;
                script += `}\n`;
            }
            script += `cat("Loaded", length(gene_sets), "gene sets from MSigDB\\n")\n\n`;
        }

        // Add custom sets if any (embedded with dictionary encoding)
        if (hasCustomSets) {
            const customGenes = new Set();
            for (const setGenes of Object.values(customSets)) {
                for (const g of setGenes) customGenes.add(g.toUpperCase());
            }
            const customDict = [...customGenes];
            const customToIdx = new Map(customDict.map((g, i) => [g, i + 1]));

            script += `# --- Custom gene sets (${Object.keys(customSets).length} sets, dictionary-encoded) ---\n`;
            script += `custom_dict <- c(${customDict.map(g => `"${this._rEscape(g)}"`).join(',')})\n`;
            script += `custom_sets <- list(\n`;
            const setEntries = [];
            for (const [name, setGenes] of Object.entries(customSets)) {
                const indices = setGenes.map(g => customToIdx.get(g.toUpperCase())).filter(i => i !== undefined);
                setEntries.push(`    "${this._rEscape(name)}" = custom_dict[c(${indices.join(',')})]`);
            }
            script += setEntries.join(',\n') + '\n)\n';
            if (useMsigdbr) {
                script += `gene_sets <- c(gene_sets, custom_sets)\n\n`;
            } else {
                script += `gene_sets <- custom_sets\n\n`;
            }
        }

        // If only custom sets and no msigdbr
        if (!useMsigdbr && !hasCustomSets) {
            // Fallback: embed all gene sets with dictionary encoding
            const allGenes = new Set(genes.map(g => g.toUpperCase()));
            for (const setGenes of Object.values(geneSets)) {
                for (const g of setGenes) allGenes.add(g.toUpperCase());
            }
            const geneDict = [...allGenes];
            const geneToIdx = new Map(geneDict.map((g, i) => [g, i + 1]));

            script += `# --- Gene sets (${Object.keys(geneSets).length} sets, dictionary-encoded) ---\n`;
            script += `gene_dict <- c(${geneDict.map(g => `"${this._rEscape(g)}"`).join(',')})\n`;
            script += `gene_sets <- list(\n`;
            const setEntries = [];
            for (const [name, setGenes] of Object.entries(geneSets)) {
                const indices = setGenes.map(g => geneToIdx.get(g.toUpperCase())).filter(i => i !== undefined);
                setEntries.push(`    "${this._rEscape(name)}" = gene_dict[c(${indices.join(',')})]`);
            }
            script += setEntries.join(',\n') + '\n)\n\n';
        }

        script += `# --- Run fgsea ---
nSets <- length(gene_sets)
cat("Running fgsea with", nSets, "gene sets and ${permutations} permutations...\\n")
results <- fgsea(
    pathways = gene_sets,
    stats = ranked_stats,
    minSize = ${minSize},
    maxSize = ${maxSize},
    nPermSimple = max(${permutations}, 10000)
)
cat("Done! Found", sum(results$padj < 0.25, na.rm = TRUE), "significant sets (FDR < 0.25)\\n")
n_na <- sum(is.na(results$padj))
if (n_na > 0) cat("Note:", n_na, "pathways had NA p-values (unbalanced gene-level stats)\\n")

# --- Format for Enrich import ---
enrich_results <- lapply(seq_len(nrow(results)), function(i) {
    row <- results[i, ]
    list(
        name = row$pathway,
        es = row$ES,
        nes = row$NES,
        pvalue = row$pval,
        fdr = row$padj,
        size = row$size,
        leadingEdge = row$leadingEdge[[1]],
        source = "fgsea"
    )
})

# --- Save as JSON ---
output_file <- "enrich_fgsea_results.json"
writeLines(toJSON(enrich_results, auto_unbox = TRUE, pretty = TRUE), output_file)
cat("Results saved to:", output_file, "\\n")
cat("Upload this file to Enrich to visualize the results.\\n")
`;

        const blob = new Blob([script], { type: 'text/plain' });
        const sizeMB = (blob.size / (1024 * 1024)).toFixed(1);
        const sizeKB = (blob.size / 1024).toFixed(0);
        const sizeStr = blob.size > 1024 * 1024 ? `${sizeMB} MB` : `${sizeKB} KB`;
        this.downloadBlob(blob, 'enrich_fgsea_analysis.R');
        this.showStatus('runStatus', 'success',
            `R script downloaded (${sizeStr}, ${genes.length} genes, ${Object.keys(geneSets).length} gene sets${useMsigdbr ? ' via msigdbr' : ''}). Run in R, then upload the JSON output here.`);
    }

    // --------------------------------------------------------
    // R Results Upload
    // --------------------------------------------------------
    async loadRResults(file) {
        if (!file) return;
        this.showStatus('uploadStatus', 'info', 'Loading fgsea results...');
        try {
            const text = await file.text();
            const data = JSON.parse(text);

            if (!Array.isArray(data) || !data[0]?.name || data[0]?.nes === undefined) {
                throw new Error('Invalid format. Expected JSON array with name, nes, pvalue, fdr fields.');
            }

            this.results = data.map(r => ({
                name: r.name,
                es: r.es || 0,
                nes: r.nes,
                pvalue: r.pvalue || r.pval || 0,
                fdr: r.fdr || r.padj || 1,
                size: r.size || 0,
                leadingEdge: r.leadingEdge || [],
                runningES: [],
                hits: [],
                source: 'fgsea'
            }));

            this.analysisDate = new Date();
            this._pinnedBubbleSets = new Set();

            const nSig = this.results.filter(r => r.fdr < 0.25).length;
            const doneMsg = `Loaded ${this.results.length} fgsea results \u2014 ${nSig} significant (FDR < 0.25)`;
            this.showStatus('uploadStatus', 'success', doneMsg);
            await this._renderResultsAsync(doneMsg);
        } catch (err) {
            this.showStatus('uploadStatus', 'error', `Failed to load R results: ${err.message}`);
        }
    }

    // --------------------------------------------------------
    // Smart Run — iterative GSEA
    // --------------------------------------------------------
    async runSmartGSEA() {
        this.buildRankedList();
        if (this.rankedList.genes.length === 0) {
            this.showStatus('runStatus', 'error', 'No valid gene-metric pairs found.');
            return;
        }
        const allGeneSets = this.getActiveGeneSets();
        if (!allGeneSets) {
            this.showStatus('runStatus', 'error', 'No gene sets loaded.');
            return;
        }
        this.readSettings();
        this._allRunGeneSets = allGeneSets;
        this._smartRunCancelled = false;

        const nTotal = Object.keys(allGeneSets).length;
        if (nTotal <= 200) {
            // Small enough to run directly
            this.runGSEA();
            return;
        }

        // Show styled dialog instead of confirm()
        const action = await this._showSmartRunDialog(nTotal);
        if (action === 'cancel') return;

        // Show progress immediately
        document.getElementById('runBtn').style.display = 'none';
        document.getElementById('smartRunBtn').style.display = 'none';
        document.getElementById('cancelBtn').style.display = '';
        document.getElementById('progressContainer').classList.add('active');
        document.getElementById('progressText').textContent = 'Initializing Smart Run...';
        document.getElementById('progressBar').style.width = '2%';
        this.hideStatus('runStatus');
        this.analysisDate = new Date();

        // Yield to let UI update before heavy computation
        await new Promise(r => setTimeout(r, 50));

        // --- Phase 1: Screen diverse sets ---
        document.getElementById('progressText').textContent = 'Phase 1/3: Computing diversity filter...';
        await new Promise(r => setTimeout(r, 0));

        const diverseSets = this._computeDiverseSubset(allGeneSets, 0.1);
        const nDiverse = Object.keys(diverseSets).length;

        if (this._smartRunCancelled) return;
        document.getElementById('progressText').textContent = `Phase 1/3: Screening ${nDiverse} diverse sets with 100 permutations...`;

        await this._runWorkerAsync(diverseSets, {
            permutations: 100,
            minSize: this.settings.minSize,
            maxSize: this.settings.maxSize,
            weightP: this.settings.weightP
        });

        if (this._smartRunCancelled) return;

        // Store Phase 1 results
        this.results = this._lastWorkerResults;
        const phase1Hits = this.results.filter(r => r.fdr < 0.25);

        // --- Phase 2: Expand around hits ---
        if (phase1Hits.length > 0) {
            document.getElementById('progressText').textContent = `Phase 2/3: Finding sets related to ${phase1Hits.length} hits...`;
            await new Promise(r => setTimeout(r, 0));

            const expansionSets = this._findRelatedSets(phase1Hits, allGeneSets, 0.1);
            const nExpansion = Object.keys(expansionSets).length;

            if (nExpansion > 0 && !this._smartRunCancelled) {
                document.getElementById('progressText').textContent = `Phase 2/3: Testing ${nExpansion} related sets...`;
                this._rerunSetNames = new Set(Object.keys(expansionSets));

                await this._runWorkerAsync(expansionSets, {
                    permutations: 100,
                    minSize: this.settings.minSize,
                    maxSize: this.settings.maxSize,
                    weightP: this.settings.weightP
                });

                if (this._smartRunCancelled) return;

                // Merge Phase 2 results
                const newResults = this._lastWorkerResults;
                const existingNames = new Set(this.results.map(r => r.name));
                for (const r of newResults) {
                    if (!existingNames.has(r.name)) this.results.push(r);
                }
                this._rerunSetNames = null;
            }
        }

        // --- Phase 3: Refine hits with high permutations ---
        const allHits = this.results.filter(r => r.fdr < 0.25);
        if (allHits.length > 0 && !this._smartRunCancelled) {
            const refinePerm = Math.max(1000, this.settings.permutations);
            document.getElementById('progressText').textContent = `Phase 3/3: Refining ${allHits.length} hits with ${refinePerm} permutations...`;

            const refineSets = {};
            for (const r of allHits) {
                if (allGeneSets[r.name]) refineSets[r.name] = allGeneSets[r.name];
            }
            this._rerunSetNames = new Set(Object.keys(refineSets));

            await this._runWorkerAsync(refineSets, {
                permutations: refinePerm,
                minSize: this.settings.minSize,
                maxSize: this.settings.maxSize,
                weightP: this.settings.weightP
            });

            if (this._smartRunCancelled) return;

            // Merge refined results
            const refinedByName = {};
            for (const r of this._lastWorkerResults) refinedByName[r.name] = r;
            this.results = this.results.map(r => refinedByName[r.name] || r);
            this._rerunSetNames = null;
        }

        // Done — render results
        document.getElementById('runBtn').style.display = '';
        document.getElementById('smartRunBtn').style.display = '';
        document.getElementById('cancelBtn').style.display = 'none';
        document.getElementById('rerunSection').style.display = '';
        this.updateRerunHitCount();

        const nSig = this.results.filter(r => r.fdr < 0.25).length;
        const nTested = this.results.length;
        const nSkipped = nTotal - nTested;
        let doneMsg = `Smart Run complete! ${nTested} gene sets tested, ${nSig} significant (FDR < 0.25)`;
        if (nSkipped > 0) {
            doneMsg += `. ${nSkipped} redundant sets skipped (similar to tested sets, Jaccard > 0.1)`;
        }
        await this._renderResultsAsync(doneMsg);
    }

    /** Run worker and return a promise that resolves when complete */
    _runWorkerAsync(geneSets, settings) {
        return new Promise((resolve, reject) => {
            const handler = (e) => {
                const data = e.data;
                if (data.type === 'progress') {
                    // Prepend phase info to progress text
                    const phasePrefix = document.getElementById('progressText').textContent.match(/^Phase \d\/\d:/)?.[0] || '';
                    document.getElementById('progressBar').style.width = data.percent + '%';
                    if (phasePrefix) {
                        document.getElementById('progressText').textContent = `${phasePrefix} ${data.text}`;
                    }
                }
                if (data.type === 'complete') {
                    this.worker.onmessage = (e) => this.handleWorkerMessage(e);
                    this._lastWorkerResults = data.results;
                    resolve(data.results);
                }
                if (data.type === 'error') {
                    this.worker.onmessage = (e) => this.handleWorkerMessage(e);
                    reject(new Error(data.message));
                }
            };
            this.worker.onmessage = handler;
            this.worker.postMessage({
                type: 'run',
                rankedGenes: this.rankedList.genes,
                rankedMetrics: this.rankedList.metrics,
                geneSets,
                settings
            });
        });
    }

    /** Compute a diverse subset of gene sets using greedy Jaccard filtering */
    _computeDiverseSubset(geneSets, maxJaccard) {
        const entries = Object.entries(geneSets);
        const kept = {};
        const keptGeneSets = [];
        for (const [name, genes] of entries) {
            const geneSet = new Set(genes.map(g => g.toUpperCase()));
            let tooSimilar = false;
            for (const keptGS of keptGeneSets) {
                const sA = geneSet.size, sB = keptGS.size;
                if (Math.min(sA, sB) / Math.max(sA, sB) < maxJaccard) continue;
                let intersection = 0;
                const smaller = sA <= sB ? geneSet : keptGS;
                const larger = sA <= sB ? keptGS : geneSet;
                for (const g of smaller) {
                    if (larger.has(g)) intersection++;
                }
                if ((sA + sB - intersection) > 0 && intersection / (sA + sB - intersection) >= maxJaccard) {
                    tooSimilar = true;
                    break;
                }
            }
            if (!tooSimilar) {
                kept[name] = genes;
                keptGeneSets.push(geneSet);
            }
        }
        return kept;
    }

    /** Find gene sets related (Jaccard > threshold) to any hit */
    _findRelatedSets(hits, allGeneSets, minJaccard) {
        const hitGeneSets = hits.map(r => ({
            name: r.name,
            genes: new Set((allGeneSets[r.name] || []).map(g => g.toUpperCase()))
        }));
        const alreadyTested = new Set(this.results.map(r => r.name));
        const related = {};
        for (const [name, genes] of Object.entries(allGeneSets)) {
            if (alreadyTested.has(name)) continue;
            const geneSet = new Set(genes.map(g => g.toUpperCase()));
            for (const hit of hitGeneSets) {
                const sA = geneSet.size, sB = hit.genes.size;
                if (Math.min(sA, sB) / Math.max(sA, sB) < minJaccard) continue;
                let intersection = 0;
                const smaller = sA <= sB ? geneSet : hit.genes;
                const larger = sA <= sB ? hit.genes : geneSet;
                for (const g of smaller) {
                    if (larger.has(g)) intersection++;
                }
                const union = sA + sB - intersection;
                if (union > 0 && intersection / union >= minJaccard) {
                    related[name] = genes;
                    break;
                }
            }
        }
        return related;
    }

    // --------------------------------------------------------
    // Expand Around Hits
    // --------------------------------------------------------
    expandAroundHits() {
        if (!this.results || !this._allRunGeneSets) {
            this.showStatus('runStatus', 'error', 'Run GSEA first to generate results.');
            return;
        }
        const hits = this.results.filter(r => r.fdr < 0.25);
        if (hits.length === 0) {
            this.showStatus('runStatus', 'warning', 'No significant hits (FDR < 0.25) to expand around.');
            return;
        }

        const related = this._findRelatedSets(hits, this._allRunGeneSets, 0.1);
        const nRelated = Object.keys(related).length;
        if (nRelated === 0) {
            this.showStatus('runStatus', 'info', 'No additional related gene sets found.');
            return;
        }

        this.readSettings();
        const perms = this.settings.permutations;
        const proceed = confirm(
            `Found ${nRelated} gene sets similar to your ${hits.length} significant hits (Jaccard > 0.1).\n\n` +
            `Run with ${perms} permutations?\n\n` +
            `This helps discover the most specific pathway driving the enrichment.`
        );
        if (!proceed) return;

        // Use the re-run mechanism to merge results
        this._rerunSetNames = new Set(Object.keys(related));
        document.getElementById('runBtn').style.display = 'none';
        document.getElementById('cancelBtn').style.display = '';
        document.getElementById('progressContainer').classList.add('active');
        document.getElementById('progressBar').style.width = '0%';
        document.getElementById('progressText').textContent = `Expanding: testing ${nRelated} related gene sets...`;

        this.worker.postMessage({
            type: 'run',
            rankedGenes: this.rankedList.genes,
            rankedMetrics: this.rankedList.metrics,
            geneSets: related,
            settings: {
                permutations: perms,
                minSize: this.settings.minSize,
                maxSize: this.settings.maxSize,
                weightP: this.settings.weightP
            }
        });
    }

    rerunFiltered() {
        if (!this.results || !this.rankedList) {
            alert('No results to re-run. Run GSEA first.');
            return;
        }

        // Get the currently filtered gene set names (same logic as filterAndRenderTable)
        const query = document.getElementById('tableFilter').value.toLowerCase();
        const fdrVal = document.getElementById('fdrFilter').value;
        const pvalVal = document.getElementById('pvalueFilter').value;
        const dirVal = document.getElementById('directionFilter').value;

        let filtered = this.results.slice();
        // Apply hidden sets filter
        if (this._hiddenSets.size > 0) {
            filtered = filtered.filter(r => !this._hiddenSets.has(r.name));
        }
        if (query) filtered = filtered.filter(r => r.name.toLowerCase().includes(query));
        if (fdrVal !== 'all') {
            const thresh = parseFloat(fdrVal);
            filtered = filtered.filter(r => r.fdr < thresh);
        }
        if (pvalVal !== 'all') {
            const thresh = parseFloat(pvalVal);
            filtered = filtered.filter(r => r.pvalue < thresh);
        }
        if (dirVal === 'up') filtered = filtered.filter(r => r.nes > 0);
        else if (dirVal === 'down') filtered = filtered.filter(r => r.nes < 0);

        if (filtered.length === 0) {
            alert('No gene sets match the current filters. Adjust filters to select gene sets for re-analysis.');
            return;
        }

        // Suggest permutations based on previous run and gene set count
        const prevPerms = this.settings.permutations;
        let suggestedPerms;
        if (prevPerms <= 500) {
            suggestedPerms = 1000;
        } else {
            suggestedPerms = prevPerms;
        }
        if (filtered.length > 100 && suggestedPerms > 500) {
            suggestedPerms = Math.min(suggestedPerms, 1000);
        }

        // Show custom dialog
        const dialog = document.getElementById('rerunDialog');
        const backdrop = document.getElementById('rerunDialogBackdrop');
        document.getElementById('rerunDialogInfo').textContent =
            `${filtered.length} gene sets match the current filters.`;
        const permInput = document.getElementById('rerunPermInput');
        permInput.value = suggestedPerms;
        document.getElementById('rerunPermHint').textContent =
            prevPerms <= 500
                ? `Previous run used ${prevPerms} permutations — suggest increasing to ≥1000 for accurate p-values.`
                : `Previous run used ${prevPerms.toLocaleString()} permutations.`;

        const updateEstimate = () => {
            const p = parseInt(permInput.value) || 1000;
            const est = Math.max(1, Math.round(filtered.length * p / 500000));
            document.getElementById('rerunEstimate').textContent =
                `Estimated time: ~${est} minute${est > 1 ? 's' : ''}`;
        };
        updateEstimate();
        permInput.oninput = updateEstimate;

        dialog.style.display = 'block';
        backdrop.style.display = 'block';

        // Wire up run button (replace handler to avoid stacking)
        const runBtn = document.getElementById('rerunDialogRunBtn');
        const newRunBtn = runBtn.cloneNode(true);
        runBtn.parentNode.replaceChild(newRunBtn, runBtn);
        newRunBtn.addEventListener('click', () => {
            dialog.style.display = 'none';
            backdrop.style.display = 'none';
            this._executeRerun(filtered, parseInt(permInput.value) || 1000);
        });
    }

    _executeRerun(filtered, permutations) {
        // Build a gene set dict from only the filtered results
        const activeGeneSets = this.getActiveGeneSets();
        const rerunSets = {};
        for (const r of filtered) {
            if (activeGeneSets[r.name]) {
                rerunSets[r.name] = activeGeneSets[r.name];
            }
        }

        if (Object.keys(rerunSets).length === 0) {
            alert('Could not find gene set data for the filtered results. The gene set collections may have changed since the last run.');
            return;
        }

        // Track which sets are being re-run so we can merge results
        this._rerunSetNames = new Set(Object.keys(rerunSets));

        // Run GSEA with only the filtered sets
        this.readSettings();
        const s = this.settings;

        document.getElementById('runBtn').style.display = 'none';
        document.getElementById('cancelBtn').style.display = '';
        document.getElementById('progressContainer').classList.add('active');
        document.getElementById('progressBar').style.width = '0%';
        document.getElementById('progressText').textContent = `Re-running ${Object.keys(rerunSets).length} gene sets...`;

        this.worker.postMessage({
            type: 'run',
            rankedGenes: this.rankedList.genes,
            rankedMetrics: this.rankedList.metrics,
            geneSets: rerunSets,
            settings: {
                permutations: permutations,
                minSize: s.minSize,
                maxSize: s.maxSize,
                weightP: s.weightP
            }
        });
    }

    rerunSignificant() {
        if (!this.results || !this.rankedList) {
            alert('No results to re-run. Run GSEA first.');
            return;
        }

        const fdrThresh = parseFloat(document.getElementById('rerunFdrThreshold').value);
        const hits = this.results.filter(r => r.fdr < fdrThresh);

        if (hits.length === 0) {
            alert(`No gene sets with FDR < ${fdrThresh}. Try a less strict threshold.`);
            return;
        }

        this.readSettings();
        const perms = this.settings.permutations;

        const proceed = confirm(
            `Re-run ${hits.length} gene sets (FDR < ${fdrThresh}) with ${perms.toLocaleString()} permutations?\n\n` +
            `This refines p-values and FDR for your hits. Use ≥1000 permutations for publication-quality statistics.\n\nContinue?`
        );
        if (!proceed) return;

        // Build gene set dict from only the significant results
        const activeGeneSets = this.getActiveGeneSets();
        const rerunSets = {};
        for (const r of hits) {
            if (activeGeneSets[r.name]) {
                rerunSets[r.name] = activeGeneSets[r.name];
            }
        }

        if (Object.keys(rerunSets).length === 0) {
            alert('Could not find gene set data for the hits. The collections may have changed since the last run.');
            return;
        }

        const s = this.settings;
        document.getElementById('runBtn').style.display = 'none';
        document.getElementById('rerunSection').style.display = 'none';
        document.getElementById('cancelBtn').style.display = '';
        document.getElementById('progressContainer').classList.add('active');
        document.getElementById('progressBar').style.width = '0%';
        document.getElementById('progressText').textContent = `Re-running ${Object.keys(rerunSets).length} gene sets...`;

        this.worker.postMessage({
            type: 'run',
            rankedGenes: this.rankedList.genes,
            rankedMetrics: this.rankedList.metrics,
            geneSets: rerunSets,
            settings: {
                permutations: s.permutations,
                minSize: s.minSize,
                maxSize: s.maxSize,
                weightP: s.weightP
            }
        });
    }

    updateRerunHitCount() {
        if (!this.results) return;
        const fdrThresh = parseFloat(document.getElementById('rerunFdrThreshold').value);
        const count = this.results.filter(r => r.fdr < fdrThresh).length;
        document.getElementById('rerunHitCount').textContent = `(${count} set${count !== 1 ? 's' : ''})`;
    }

    handleWorkerMessage(e) {
        const data = e.data;

        if (data.type === 'progress') {
            document.getElementById('progressBar').style.width = data.percent + '%';
            document.getElementById('progressText').textContent = data.text;
            // Capture size filter info for completion message
            if (data.sizeFilterInfo) this._lastSizeFilterInfo = data.sizeFilterInfo;
        }

        if (data.type === 'complete') {
            // If this was a re-run/expand, merge new results into existing results
            if (this._rerunSetNames && this._rerunSetNames.size > 0 && this.results) {
                const newResultsByName = {};
                for (const r of data.results) newResultsByName[r.name] = r;
                // Update existing results
                const existingNames = new Set();
                this.results = this.results.map(r => {
                    existingNames.add(r.name);
                    return newResultsByName[r.name] ? newResultsByName[r.name] : r;
                });
                // Add any new results not previously in the list
                for (const r of data.results) {
                    if (!existingNames.has(r.name)) this.results.push(r);
                }
                this._rerunSetNames = null;
            } else {
                this.results = data.results;
            }
            document.getElementById('runBtn').style.display = '';
            document.getElementById('smartRunBtn').style.display = '';
            document.getElementById('cancelBtn').style.display = 'none';

            const nSig = this.results.filter(r => r.fdr < 0.25).length;
            let doneMsg = `Done! ${this.results.length} gene sets tested, ${nSig} significant (FDR < 0.25)`;
            // Add exclusion details if any sets were filtered out
            const sfi = this._lastSizeFilterInfo;
            if (sfi && sfi.totalInput > sfi.totalPassed) {
                const excluded = sfi.totalInput - sfi.totalPassed;
                const reasons = [];
                if (sfi.excludedTooSmall > 0) reasons.push(`${sfi.excludedTooSmall} had fewer than ${this.settings.minSize} matching genes (Min size)`);
                if (sfi.excludedTooLarge > 0) reasons.push(`${sfi.excludedTooLarge} had more than ${this.settings.maxSize} matching genes (Max size)`);
                if (sfi.excludedNoOverlap > 0) reasons.push(`${sfi.excludedNoOverlap} had no genes in your data`);
                doneMsg += `. ${excluded} excluded (${reasons.join(', ')})`;
            }
            this._lastSizeFilterInfo = null;

            // Show re-run options now that we have results
            document.getElementById('rerunSection').style.display = '';
            this.updateRerunHitCount();

            // Render results asynchronously with step-by-step progress
            this._renderResultsAsync(doneMsg);
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
    // Async Rendering Pipeline
    // --------------------------------------------------------
    async _renderResultsAsync(doneMsg) {
        const progressBar = document.getElementById('progressBar');
        const progressText = document.getElementById('progressText');
        progressBar.style.width = '100%';

        // Elapsed timer
        const startTime = Date.now();
        const timer = setInterval(() => {
            const elapsed = Math.round((Date.now() - startTime) / 1000);
            const base = this._renderStepLabel || 'Rendering results';
            progressText.textContent = `${base}... (${elapsed}s)`;
        }, 1000);

        const step = async (label, fn) => {
            this._renderStepLabel = label;
            progressText.textContent = `${label}...`;
            await new Promise(r => setTimeout(r, 0)); // yield to browser
            try { fn(); } catch (err) { console.error(`Error in ${label}:`, err); }
        };

        try {
            // Show result panels
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
            this._syncSidebarOffset();

            await step('Rendering overview', () => {
                this.populateGeneSetSelector();
                this.renderBubblePlot();
            });
            await step('Rendering ranked list', () => {
                this.renderRankedPlot();
            });
            await step('Building results table', () => {
                this.filterAndRenderTable();
                this.generateMethods();
                this.renderTopBottomGenes();
            });
            await step('Rendering enrichment plot', () => {
                if (this.results.length > 0) {
                    const topSet = this.results[0].name;
                    document.getElementById('geneSetSelector').value = topSet;
                    this.renderESPlot(topSet);
                    this.renderGeneSetInfo(topSet);
                    this.renderGeneDetailTable(topSet);
                }
            });
            await step('Building overlap data', () => {
                this._buildOverlapCache();
                this.renderOverlapHeatmap();
            });
        } finally {
            clearInterval(timer);
            this._renderStepLabel = null;
            document.getElementById('progressContainer').classList.remove('active');
            this.showStatus('runStatus', 'success', doneMsg);
        }
    }

    // --------------------------------------------------------
    // Display Results (called for re-renders, not initial load)
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

        // Align sidebar top with main panel (account for tab bar height)
        this._syncSidebarOffset();

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

    populateGeneSetSelector(searchFilter) {
        const sel = document.getElementById('geneSetSelector');
        const currentVal = sel.value;
        sel.innerHTML = '<option value="">Select a gene set...</option>';

        // Read enrichment plot tab filters
        const esFdrEl = document.getElementById('esFdrFilter');
        const esPvalEl = document.getElementById('esPvalFilter');
        const esDirEl = document.getElementById('esDirectionFilter');
        const esRedEl = document.getElementById('esRedundancyFilter');
        const fdrThresh = esFdrEl ? parseFloat(esFdrEl.value) : 1;
        const pvalThresh = esPvalEl ? parseFloat(esPvalEl.value) : 1;
        const dirFilter = esDirEl ? esDirEl.value : 'all';
        const redundancyThresh = esRedEl ? parseFloat(esRedEl.value) : 0;

        // Filter results by FDR/pval/direction + hidden sets
        let filtered = this.results.filter(r => {
            if (this._hiddenSets.has(r.name)) return false;
            if (fdrThresh < 1 && r.fdr >= fdrThresh) return false;
            if (pvalThresh < 1 && r.pvalue >= pvalThresh) return false;
            if (dirFilter === 'up' && r.nes <= 0) return false;
            if (dirFilter === 'down' && r.nes > 0) return false;
            return true;
        });

        // Apply redundancy filter — keep only cluster representatives
        if (redundancyThresh > 0 && filtered.length > 1) {
            const clusters = this._computeOverlapClusters(filtered, redundancyThresh);
            const reps = this._getClusterRepresentatives(clusters);
            // Keep sets that are representatives OR not in any cluster
            filtered = filtered.filter(r => reps.has(r.name) || !clusters[r.name]);
        }

        // Capture count before search filter
        const totalFiltered = filtered.length;

        // Apply search filter if provided
        if (searchFilter) {
            const q = searchFilter.toLowerCase();
            filtered = filtered.filter(r => r.name.toLowerCase().includes(q));
        }

        // Update count
        const countEl = document.getElementById('esFilterCount');
        if (countEl) {
            countEl.textContent = searchFilter
                ? `${filtered.length} matches (${totalFiltered} of ${this.results.length} gene sets)`
                : `${totalFiltered} of ${this.results.length} gene sets`;
        }

        // Group by positive/negative NES
        const pos = filtered.filter(r => r.nes > 0).sort((a, b) => b.nes - a.nes);
        const neg = filtered.filter(r => r.nes <= 0).sort((a, b) => a.nes - b.nes);

        if (pos.length > 0) {
            const group = document.createElement('optgroup');
            group.label = `▲ Upregulated (${pos.length})`;
            for (const r of pos) {
                const opt = new Option(
                    `${this.cleanName(r.name)}  —  NES: ${r.nes.toFixed(2)}, FDR: ${this.formatPval(r.fdr)}`,
                    r.name
                );
                group.appendChild(opt);
            }
            sel.appendChild(group);
        }

        if (neg.length > 0) {
            const group = document.createElement('optgroup');
            group.label = `▼ Downregulated (${neg.length})`;
            for (const r of neg) {
                const opt = new Option(
                    `${this.cleanName(r.name)}  —  NES: ${r.nes.toFixed(2)}, FDR: ${this.formatPval(r.fdr)}`,
                    r.name
                );
                group.appendChild(opt);
            }
            sel.appendChild(group);
        }

        // Restore previous selection if still exists
        if (currentVal) {
            const stillExists = Array.from(sel.options).some(o => o.value === currentVal);
            if (stillExists) sel.value = currentVal;
        }
    }

    viewPinnedInOverview() {
        if (this._pinnedBubbleSets.size === 0) return;
        this._showOnlyPinned = true;
        this.renderBubblePlot();
        this.showTab('overview');
    }

    // --------------------------------------------------------
    // Bubble Plot
    // --------------------------------------------------------
    renderBubblePlot() {
        this.readSettings();
        // Read from overview toolbar (primary) with gear popup as fallback
        const ovFdr = document.getElementById('overviewFdrFilter');
        const ovPval = document.getElementById('overviewPvalFilter');
        const ovTopN = document.getElementById('overviewTopN');
        const ovDir = document.getElementById('overviewDirectionFilter');
        const ovOverlap = document.getElementById('overviewOverlapFilter');
        const fdrThresh = ovFdr ? parseFloat(ovFdr.value) : parseFloat(this.settings.fdrDisplayThreshold);
        const pvalThresh = ovPval ? parseFloat(ovPval.value) : parseFloat(document.getElementById('pvalDisplayThreshold')?.value || '1');
        const topN = ovTopN ? parseInt(ovTopN.value) : this.settings.topN;
        const dirFilter = ovDir ? ovDir.value : 'all';
        const overlapThresh = ovOverlap ? parseInt(ovOverlap.value) : this._overlapFilterThreshold;
        const fontFam = this.settings.fontFamily + ', sans-serif';
        const baseFontSize = this.settings.fontSize;

        let filtered = this.results;
        // Apply hidden sets filter
        if (this._hiddenSets.size > 0) {
            filtered = filtered.filter(r => !this._hiddenSets.has(r.name));
        }
        if (fdrThresh < 1) {
            filtered = filtered.filter(r => r.fdr < fdrThresh);
        }
        if (pvalThresh < 1) {
            filtered = filtered.filter(r => r.pvalue < pvalThresh);
        }
        // Direction filter
        if (dirFilter === 'up') {
            filtered = filtered.filter(r => r.nes > 0);
        } else if (dirFilter === 'down') {
            filtered = filtered.filter(r => r.nes < 0);
        }
        // Apply overlap redundancy filter
        if (overlapThresh > 0) {
            filtered = this._applyOverlapFilter(filtered, overlapThresh);
        }
        // Update result count
        const countEl = document.getElementById('overviewResultCount');
        if (countEl) countEl.textContent = `${filtered.length} sets`;
        // Sort by |NES| and take top N, plus any pinned sets
        const sorted = filtered.sort((a, b) => Math.abs(b.nes) - Math.abs(a.nes));
        let top;
        if (this._showOnlyPinned && this._pinnedBubbleSets.size > 0) {
            // Show ONLY pinned sets (bypass all filters)
            top = this.results.filter(r => this._pinnedBubbleSets.has(r.name));
            top.sort((a, b) => Math.abs(b.nes) - Math.abs(a.nes));
            this._showOnlyPinned = false; // reset flag after use
        } else {
            const topSet = new Set(sorted.slice(0, topN).map(r => r.name));
            // Add pinned sets — bypass all filters (use unfiltered results)
            for (const r of this.results) {
                if (this._pinnedBubbleSets.has(r.name) && !topSet.has(r.name)) {
                    topSet.add(r.name);
                }
            }
            top = sorted.filter(r => topSet.has(r.name));
            // Also add pinned sets that were filtered out
            for (const r of this.results) {
                if (this._pinnedBubbleSets.has(r.name) && !top.find(t => t.name === r.name)) {
                    top.push(r);
                }
            }
        }

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
        // Plot border lines at paper coordinates (top + bottom only — no left/right for a clean open look)
        shapes.push(
            { type: 'line', xref: 'paper', yref: 'paper', x0: 0, x1: 1, y0: 1, y1: 1, line: { color: '#333', width: 1 } },  // top
            { type: 'line', xref: 'paper', yref: 'paper', x0: 0, x1: 1, y0: 0, y1: 0, line: { color: '#333', width: 1 } }   // bottom
        );

        const cbTitleFont = this._getTextFont('bubble', 'colorbarTitle');
        const cbTickFont = this._getTextFont('bubble', 'colorbarTickFont');
        const trace = {
            x: top.map(r => r.nes),
            y: top.map((_, i) => i),
            mode: 'markers',
            type: 'scatter',
            marker: {
                size: top.map(r => Math.max(10, Math.min(30, Math.sqrt(r.size) * 2))),
                color: fdrColors,
                colorscale: this.settings.colorScale,
                cmin: 0,
                cmax: maxLogFdr,
                showscale: true,
                colorbar: {
                    title: { text: cbTitleFont.visible !== false ? cbTitleFont.wrap(cbTitleFont.text || 'FDR') : '', font: { size: cbTitleFont.size, family: cbTitleFont.family } },
                    thickness: 20,
                    len: 0.35,
                    y: 0.8,
                    x: 1.02,
                    tickvals: [0, -Math.log10(0.25), -Math.log10(0.1), -Math.log10(0.05), -Math.log10(0.01)].filter(v => v <= maxLogFdr),
                    ticktext: ['1', '0.25', '0.1', '0.05', '0.01'].slice(0, [0, -Math.log10(0.25), -Math.log10(0.1), -Math.log10(0.05), -Math.log10(0.01)].filter(v => v <= maxLogFdr).length),
                    tickfont: { size: cbTickFont.visible !== false ? cbTickFont.size : 0, family: cbTickFont.family || fontFam },
                    outlinewidth: 1,
                    outlinecolor: '#ccc'
                },
                line: { width: 1, color: 'rgba(0,0,0,0.3)' }
            },
            text: top.map(r =>
                `<b>${this.cleanName(r.name)}</b><br>` +
                `NES: ${r.nes.toFixed(3)}<br>` +
                `FDR: ${this.formatPval(r.fdr)}<br>` +
                `p-value: ${this.formatPval(r.pvalue)}<br>` +
                `Size: ${r.size} genes`
            ),
            hoverinfo: 'text',
            cliponaxis: false,
            showlegend: false
        };

        const bubbleXFont = this._getTextFont('bubble', 'xAxisLabel');
        const bubbleYTickFont = this._getTextFont('bubble', 'yTickFont');
        const bubbleTitleFont = this._getTextFont('bubble', 'title');
        const layout = {
            xaxis: {
                zeroline: false,
                gridcolor: '#f0f0f0',
                side: 'bottom',
                tickfont: { size: bubbleYTickFont.visible !== false ? bubbleYTickFont.size : 0, family: bubbleYTickFont.family },
                showline: false,
                fixedrange: true
            },
            yaxis: {
                tickvals: top.map((_, i) => i),
                ticktext: top.map(r => this.cleanName(r.name)),
                tickfont: { size: bubbleYTickFont.visible !== false ? bubbleYTickFont.size : 0, family: bubbleYTickFont.family },
                automargin: true,
                gridwidth: 0,
                showgrid: false,
                zeroline: false,
                showline: false,
                range: [-1.5, top.length + 0.5],
                fixedrange: true
            },
            height: Math.max(440, top.length * 26 + 160),
            margin: { l: 15, r: 160, t: 70, b: 55 },
            font: { family: fontFam },
            paper_bgcolor: this.settings.transparentBg ? 'rgba(0,0,0,0)' : '#fff',
            plot_bgcolor: '#fff',
            shapes: shapes,
            dragmode: false
        };

        // Size annotation: Plotly caps legend symbol sizes at ~16px, so bubble legends
        // are inaccurate for large markers. Use a text annotation instead.
        const actualSizes = top.map(r => r.size).sort((a, b) => a - b);
        const minActualSize = actualSizes[0];
        const maxActualSize = actualSizes[actualSizes.length - 1];
        const sizeRangeText = minActualSize === maxActualSize
            ? `(${minActualSize} genes)`
            : `(${minActualSize} – ${maxActualSize})`;
        if (!layout.annotations) layout.annotations = [];
        // Title annotation (draggable)
        if (bubbleTitleFont.visible !== false) {
            layout.annotations.push({
                text: bubbleTitleFont.wrap(bubbleTitleFont.text || 'Enrichment Overview'),
                xref: 'paper', yref: 'paper', x: 0.5, y: 1.08,
                showarrow: false,
                font: { size: bubbleTitleFont.size, family: bubbleTitleFont.family },
                xanchor: 'center', yanchor: 'bottom'
            });
        }
        // X-axis label as annotation (draggable)
        if (bubbleXFont.visible !== false) {
            layout.annotations.push({
                text: bubbleXFont.wrap(bubbleXFont.text || 'Normalized Enrichment Score (NES)'),
                xref: 'paper', yref: 'paper', x: 0.5, y: -0.06,
                showarrow: false,
                font: { size: bubbleXFont.size, family: bubbleXFont.family },
                xanchor: 'center', yanchor: 'top'
            });
        }
        const sizeAnnotFont = this._getTextFont('bubble', 'sizeAnnotation');
        layout.annotations.push({
            text: `<b>Bubble size</b><br>= gene set size<br>${sizeRangeText}`,
            xref: 'paper', yref: 'paper',
            x: 1.02, y: 0.45,
            xanchor: 'left', yanchor: 'top',
            showarrow: false,
            font: { size: sizeAnnotFont.size, family: sizeAnnotFont.family, color: '#555' },
            align: 'center',
            visible: sizeAnnotFont.visible !== false
        });

        // Apply custom dimensions
        const dims = this.getPlotDimensions('bubblePlot', undefined, layout.height);
        if (dims.width) layout.width = dims.width;
        if (dims.height) layout.height = dims.height;

        Plotly.newPlot('bubblePlot', [trace], layout, {
            responsive: true,
            displaylogo: false,
            modeBarButtonsToRemove: ['lasso2d', 'select2d', 'zoom2d', 'pan2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'resetScale2d'],
            scrollZoom: false,
            doubleClick: false,
            edits: { annotationPosition: true, annotationText: true, axisTitleText: false, titleText: false, legendPosition: true, colorbarPosition: true }
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

        // Parse positive/negative base colors for gradient
        const posRGB = this._parseHex(this.settings.positiveColor);
        const negRGB = this._parseHex(this.settings.negativeColor);

        for (let i = 0; i < N; i += step) {
            xVals.push(i);
            yVals.push(metrics[i]);
            // Gradient: interpolate from light to full color based on magnitude
            if (metrics[i] >= 0) {
                const intensity = Math.min(1, metrics[i] / (Math.abs(metrics[0]) || 1));
                const r = Math.round(255 - (255 - posRGB[0]) * intensity);
                const g = Math.round(255 - (255 - posRGB[1]) * intensity);
                const b = Math.round(255 - (255 - posRGB[2]) * intensity);
                barColors.push(`rgba(${r}, ${g}, ${b}, 0.85)`);
            } else {
                const intensity = Math.min(1, Math.abs(metrics[i]) / (Math.abs(metrics[N - 1]) || 1));
                const r = Math.round(255 - (255 - negRGB[0]) * intensity);
                const g = Math.round(255 - (255 - negRGB[1]) * intensity);
                const b = Math.round(255 - (255 - negRGB[2]) * intensity);
                barColors.push(`rgba(${r}, ${g}, ${b}, 0.85)`);
            }
        }

        const trace = {
            x: xVals,
            y: yVals,
            type: 'bar',
            width: step * 1.05,
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
        // Title annotation (draggable)
        const rpTitleFont = this._getTextFont('ranked', 'title');
        if (rpTitleFont.visible !== false) {
            annotations.push({
                text: rpTitleFont.wrap(rpTitleFont.text || 'Ranked List Metric'),
                xref: 'paper', yref: 'paper', x: 0.5, y: 1.12,
                showarrow: false,
                font: { size: rpTitleFont.size, family: rpTitleFont.family },
                xanchor: 'center', yanchor: 'bottom'
            });
        }
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
        // Positive / Negative labels (data-type-aware)
        const rpDataType = this.settings.dataType || 'expression';
        const rpPosLabel = rpDataType === 'crispr' ? 'Enriched in screen' : 'Positively correlated';
        const rpNegLabel = rpDataType === 'crispr' ? 'Depleted in screen' : 'Negatively correlated';
        const rpPosFont = this._getTextFont('ranked', 'posLabel');
        const rpNegFont = this._getTextFont('ranked', 'negLabel');
        if (rpPosFont.visible !== false) {
            annotations.push({
                text: rpPosFont.wrap(rpPosFont.text || rpPosLabel),
                xref: 'paper', yref: 'paper', x: 0.02, y: 1.02,
                showarrow: false, font: { size: rpPosFont.size, family: rpPosFont.family, color: this.settings.positiveColor },
                xanchor: 'left', yanchor: 'bottom'
            });
        }
        if (rpNegFont.visible !== false) {
            annotations.push({
                text: rpNegFont.wrap(rpNegFont.text || rpNegLabel),
                xref: 'paper', yref: 'paper', x: 0.98, y: 1.02,
                showarrow: false, font: { size: rpNegFont.size, family: rpNegFont.family, color: this.settings.negativeColor },
                xanchor: 'right', yanchor: 'bottom'
            });
        }

        // Add range padding so extreme bars aren't clipped at edges
        const yPad = Math.max(Math.abs(metrics[0]), Math.abs(metrics[N - 1])) * 0.06;
        const xPad = N * 0.015;

        const rpXFont = this._getTextFont('ranked', 'xAxisLabel');
        const rpYFont = this._getTextFont('ranked', 'yAxisLabel');
        const rpXTickFont = this._getTextFont('ranked', 'xTickFont');
        const rpYTickFont = this._getTextFont('ranked', 'yTickFont');
        const layout = {
            xaxis: {
                title: { text: rpXFont.visible !== false ? rpXFont.wrap(rpXFont.text || 'Rank in Ordered Dataset') : '', font: { size: rpXFont.size, family: rpXFont.family } },
                showgrid: false,
                tickfont: { size: rpXTickFont.visible !== false ? rpXTickFont.size : 0, family: rpXTickFont.family },
                range: [-xPad, N - 1 + xPad],
                fixedrange: true
            },
            yaxis: {
                title: { text: rpYFont.visible !== false ? rpYFont.wrap(rpYFont.text || ('Ranked list metric (' + metricLabel + ')')) : '', font: { size: rpYFont.size, family: rpYFont.family } },
                zeroline: true,
                zerolinewidth: 1.5,
                zerolinecolor: '#333',
                gridcolor: '#e5e5e5',
                tickfont: { size: rpYTickFont.visible !== false ? rpYTickFont.size : 0, family: rpYTickFont.family },
                range: [Math.min(0, metrics[N - 1]) - yPad, Math.max(0, metrics[0]) + yPad],
                fixedrange: true
            },
            height: 250,
            margin: { l: 65, r: 20, t: 50, b: 50 },
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
            annotations: annotations,
            dragmode: false
        };

        // Apply custom dimensions
        const dims = this.getPlotDimensions('rankedPlot', undefined, layout.height);
        if (dims.width) layout.width = dims.width;
        if (dims.height) layout.height = dims.height;

        Plotly.newPlot('rankedPlot', [trace], layout, {
            responsive: true,
            displaylogo: false,
            modeBarButtonsToRemove: ['lasso2d', 'select2d', 'zoom2d', 'pan2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'resetScale2d'],
            scrollZoom: false,
            doubleClick: false,
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

        // R-imported results don't have running ES data
        if (result.source === 'fgsea' && (!result.runningES || result.runningES.length === 0)) {
            const container = document.getElementById('esPlot');
            container.innerHTML = `<div style="text-align:center; padding: 40px 20px; color: var(--gray-500); font-size: 0.9em;">
                <div style="font-weight: 500; margin-bottom: 6px;">ES curve not available for R-imported results</div>
                <div style="font-size: 0.85em;">To see the enrichment plot, re-run this gene set in the browser using "Re-run filtered".</div>
            </div>`;
            return;
        }

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

        // ---- Panel domain proportions (from settings) ----
        const pES = s.esPanelES || 55;
        const pHits = s.esPanelHits || 8;
        const pMetric = s.esPanelMetric || 28;
        const pTotal = pES + pHits + pMetric;
        const gapPct = Math.max(0, (100 - pTotal) / 2);  // Two gaps between three panels
        const usable = 0.97;  // 0.0 to 0.97 (top 3% for title)
        const scale = usable / 100;

        const metBot = 0.0;
        const metTop = pMetric * scale;
        const hitBot = (pMetric + gapPct) * scale;
        const hitTop = (pMetric + gapPct + pHits) * scale;
        const esBot = (pMetric + gapPct + pHits + gapPct) * scale;
        const esTop = usable;
        const xLeft = 0.0;
        const xRight = 1.0;

        // Range padding so extreme ends aren't clipped at borders
        const xPad = N * 0.025;

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
        // Use coarser subsampling (~600 bars) to match the ranked plot visual appearance
        const metMaxBars = 600;
        const metStep = Math.max(1, Math.ceil(N / metMaxBars));
        const metBarX = [], metBarY = [];
        for (let i = 0; i < N; i += metStep) {
            metBarX.push(i);
            metBarY.push(metrics[i]);
        }
        if (metBarX[metBarX.length - 1] !== N - 1) {
            metBarX.push(N - 1);
            metBarY.push(metrics[N - 1]);
        }

        // Build per-bar gradient colors matching the ranked plot style
        const metBarColors = [];
        if (s.showMetricFill) {
            const posRGB = this._parseHex(s.positiveColor);
            const negRGB = this._parseHex(s.negativeColor);
            const maxPos = Math.abs(metrics[0]) || 1;
            const maxNeg = Math.abs(metrics[N - 1]) || 1;
            for (let i = 0; i < metBarY.length; i++) {
                if (metBarY[i] >= 0) {
                    const intensity = Math.min(1, metBarY[i] / maxPos);
                    const r = Math.round(255 - (255 - posRGB[0]) * intensity);
                    const g = Math.round(255 - (255 - posRGB[1]) * intensity);
                    const b = Math.round(255 - (255 - posRGB[2]) * intensity);
                    metBarColors.push(`rgba(${r}, ${g}, ${b}, 0.85)`);
                } else {
                    const intensity = Math.min(1, Math.abs(metBarY[i]) / maxNeg);
                    const r = Math.round(255 - (255 - negRGB[0]) * intensity);
                    const g = Math.round(255 - (255 - negRGB[1]) * intensity);
                    const b = Math.round(255 - (255 - negRGB[2]) * intensity);
                    metBarColors.push(`rgba(${r}, ${g}, ${b}, 0.85)`);
                }
            }
        }

        const traces3 = [];
        if (s.showMetricFill) {
            // Single bar trace with per-bar gradient colors (matches ranked plot exactly)
            // Explicit width ensures bars don't collapse to 0px due to subpixel rounding
            traces3.push({
                x: metBarX, y: metBarY,
                type: 'bar',
                width: metStep * 1.05,
                marker: { color: metBarColors, line: { width: 0 } },
                showlegend: false, xaxis: 'x3', yaxis: 'y3', hoverinfo: 'skip'
            });
        } else {
            traces3.push({
                x: metBarX, y: metBarY,
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

        // ES peak indicator — vertical dashed line at the rank where ES peaks
        if (s.showESIndicator !== false) {
            // Find the position of maximum |ES| in runningES
            let esPeakIdx = 0;
            let esPeakVal = 0;
            for (let i = 0; i < result.runningES.length; i++) {
                if (Math.abs(result.runningES[i]) > Math.abs(esPeakVal)) {
                    esPeakVal = result.runningES[i];
                    esPeakIdx = i;
                }
            }
            // Vertical line spanning all three panels
            panelShapes.push({
                type: 'line', x0: esPeakIdx, x1: esPeakIdx,
                y0: 0, y1: 1, yref: 'paper',
                xref: 'x', // same x-axis as ES panel
                line: { color: '#e67e22', width: 1.2, dash: 'dash' }
            });
        }

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
        const esTitleFont = this._getTextFont('es', 'title');
        const annotations = [
            // Title
            {
                text: esTitleFont.wrap(esTitleFont.text || geneSetName.replace(/_/g, ' ')),
                xref: 'paper', yref: 'paper', x: 0.5, y: 1.06,
                showarrow: false,
                font: { size: esTitleFont.size, family: esTitleFont.family },
                xanchor: 'center',
                visible: esTitleFont.visible !== false
            }
        ];

        // Stats annotation — position in corner with least overlap with ES curve fill
        // Uses cache to prevent repositioning when settings change
        if (s.showStatsBox) {
            let bestCorner = this._statsBoxCornerCache[geneSetName];

            if (!bestCorner) {
                // First render for this gene set — compute optimal corner
                // Score each corner by how much the ES curve occupies that region
                const esValues = result.runningES;
                const esMax = Math.max(...esValues);
                const esMin = Math.min(...esValues);
                const esRange = Math.max(esMax - esMin, 0.1);

                const corners = [
                    { x: 0.03, y: esTop - (esTop - esBot) * 0.06, xa: 'left', ya: 'top' },
                    { x: 0.97, y: esTop - (esTop - esBot) * 0.06, xa: 'right', ya: 'top' },
                    { x: 0.03, y: esBot + (esTop - esBot) * 0.06, xa: 'left', ya: 'bottom' },
                    { x: 0.97, y: esBot + (esTop - esBot) * 0.06, xa: 'right', ya: 'bottom' }
                ];

                bestCorner = corners[0];
                let bestScore = Infinity;
                const boxFraction = 0.35;
                // The box occupies roughly the top or bottom 30% of the y-range
                const boxYFraction = 0.30;

                for (const corner of corners) {
                    const isTop = corner.ya === 'top';
                    const xStart = corner.xa === 'left' ? 0 : Math.floor(N * (1 - boxFraction));
                    const xEnd = corner.xa === 'left' ? Math.floor(N * boxFraction) : N - 1;
                    const sampleStep = Math.max(1, Math.floor((xEnd - xStart) / 50));

                    let overlapSum = 0;
                    let count = 0;
                    for (let i = xStart; i <= xEnd; i += sampleStep) {
                        const v = esValues[i] || 0;
                        // Normalize v to 0..1 where 1 = esMax, 0 = esMin
                        const norm = (v - esMin) / esRange;
                        // Top box occupies norm > (1 - boxYFraction), bottom occupies norm < boxYFraction
                        // Score = how much the curve is inside the box area
                        if (isTop && norm > (1 - boxYFraction)) {
                            overlapSum += (norm - (1 - boxYFraction)) / boxYFraction;
                        } else if (!isTop && norm < boxYFraction) {
                            overlapSum += (boxYFraction - norm) / boxYFraction;
                        }
                        count++;
                    }
                    const score = count > 0 ? overlapSum / count : 0;
                    if (score < bestScore) { bestScore = score; bestCorner = corner; }
                }

                this._statsBoxCornerCache[geneSetName] = bestCorner;
            } else {
                // Cached — update y coordinates for current panel proportions
                bestCorner = { ...bestCorner };
                bestCorner.y = bestCorner.ya === 'top'
                    ? esTop - (esTop - esBot) * 0.06
                    : esBot + (esTop - esBot) * 0.06;
            }

            const esStatsFont = this._getTextFont('es', 'statsBox');
            annotations.push({
                text: `NES = ${result.nes.toFixed(2)}<br>FDR = ${this.formatPval(result.fdr)}<br>p = ${this.formatPval(result.pvalue)}<br>Size = ${result.size}`,
                xref: 'paper', yref: 'paper',
                x: bestCorner.x, y: bestCorner.y,
                showarrow: false,
                font: { size: esStatsFont.size, family: esStatsFont.family || 'Roboto Mono, monospace', color: '#333' },
                align: bestCorner.xa === 'right' ? 'right' : 'left',
                xanchor: bestCorner.xa,
                yanchor: bestCorner.ya,
                bgcolor: 'rgba(255,255,255,0.95)',
                bordercolor: '#ccc',
                borderwidth: 1,
                borderpad: 5,
                visible: esStatsFont.visible !== false
            });
        }

        // Compute y-range for metric panel — use actual min/max of ALL metrics
        // (must be computed before zero-cross annotation which uses metYMin)
        let metActualMin = 0, metActualMax = 0;
        for (let i = 0; i < N; i++) {
            if (metrics[i] < metActualMin) metActualMin = metrics[i];
            if (metrics[i] > metActualMax) metActualMax = metrics[i];
        }
        const metYMin = metActualMin * 1.05;
        const metYMax = metActualMax * 1.05;

        // Zero cross annotation — positioned below the metric panel to avoid overlapping bars
        if (zeroCross >= 0 && s.showZeroCross) {
            // Place the annotation below the zero line with arrow pointing up
            annotations.push({
                text: `Zero cross at ${zeroCross.toLocaleString()}`,
                x: zeroCross, y: metYMin * 0.85,
                xref: 'x3', yref: 'y3',
                showarrow: true, arrowhead: 0, arrowsize: 0.8, arrowwidth: 1,
                ax: 0, ay: 18,
                arrowcolor: '#999',
                font: { size: Math.max(8, baseFontSize - 4), color: '#666', family: fontFam }
            });
        }

        // Correlation labels — position below the x-axis title (data-type-aware)
        if (s.showCorrelationLabels && zeroCross >= 0) {
            const esDT = this.settings.dataType || 'expression';
            const esPosLabel = esDT === 'crispr' ? 'Enriched in screen' : 'Positively correlated';
            const esNegLabel = esDT === 'crispr' ? 'Depleted in screen' : 'Negatively correlated';
            const esPosFont = this._getTextFont('es', 'posLabel');
            const esNegFont = this._getTextFont('es', 'negLabel');
            if (esPosFont.visible) {
                annotations.push({
                    text: esPosFont.wrap(esPosFont.text || esPosLabel),
                    xref: 'paper', yref: 'paper', x: 0.0, y: -0.14,
                    showarrow: false,
                    font: { size: esPosFont.size, family: esPosFont.family, color: s.positiveColor },
                    xanchor: 'left', yanchor: 'top'
                });
            }
            if (esNegFont.visible) {
                annotations.push({
                    text: esNegFont.wrap(esNegFont.text || esNegLabel),
                    xref: 'paper', yref: 'paper', x: 1.0, y: -0.14,
                    showarrow: false,
                    font: { size: esNegFont.size, family: esNegFont.family, color: s.negativeColor },
                    xanchor: 'right', yanchor: 'top'
                });
            }
        }

        // Text font references for ES plot
        const esYLabelFont = this._getTextFont('es', 'esYLabel');
        const esXFont = this._getTextFont('es', 'xAxisLabel');
        const esMetricFont = this._getTextFont('es', 'metricYLabel');
        const esTickFont = this._getTextFont('es', 'esTickFont');
        const esMetTickFont = this._getTextFont('es', 'tickFont');
        const esHitFont = this._getTextFont('es', 'hitMarkersLabel');

        // Add axis labels as draggable annotations instead of axis titles
        // Use different x positions to avoid overlap between panels
        const yLabelSize = Math.min(esYLabelFont.size, 12);  // cap size for consistency
        if (esYLabelFont.visible !== false) {
            annotations.push({
                text: esYLabelFont.wrap(esYLabelFont.text || 'Enrichment score (ES)'),
                xref: 'paper', yref: 'paper',
                x: -0.09, y: (esBot + esTop) / 2,
                showarrow: false, textangle: -90,
                font: { size: esYLabelFont.size, family: esYLabelFont.family },
                xanchor: 'center', yanchor: 'middle'
            });
        }
        if (esHitFont.visible !== false) {
            annotations.push({
                text: esHitFont.wrap(esHitFont.text || 'Gene hits'),
                xref: 'paper', yref: 'paper',
                x: -0.03, y: (hitBot + hitTop) / 2,
                showarrow: false, textangle: -90,
                font: { size: Math.min(esHitFont.size, 9), family: esHitFont.family },
                xanchor: 'center', yanchor: 'middle'
            });
        }
        if (esMetricFont.visible !== false) {
            annotations.push({
                text: esMetricFont.wrap(esMetricFont.text || metricLabel),
                xref: 'paper', yref: 'paper',
                x: -0.09, y: (metBot + metTop) / 2,
                showarrow: false, textangle: -90,
                font: { size: esMetricFont.size, family: esMetricFont.family },
                xanchor: 'center', yanchor: 'middle'
            });
        }
        if (esXFont.visible !== false) {
            annotations.push({
                text: esXFont.wrap(esXFont.text || 'Rank in Ordered Dataset'),
                xref: 'paper', yref: 'paper',
                x: 0.5, y: -0.04,
                showarrow: false,
                font: { size: esXFont.size, family: esXFont.family },
                xanchor: 'center', yanchor: 'top'
            });
        }

        const layout = {
            // ES panel
            xaxis: {
                range: [-xPad, N - 1 + xPad], showticklabels: false, showgrid: false, zeroline: false,
                domain: [xLeft, xRight], anchor: 'y', showline: false, fixedrange: true
            },
            yaxis: {
                domain: [esBot, esTop], anchor: 'x',
                gridcolor: '#eee', gridwidth: 1,
                zeroline: false,
                tickfont: { size: esTickFont.visible !== false ? esTickFont.size : 0, family: esTickFont.family || fontFam },
                fixedrange: true
            },
            // Hit marker panel
            xaxis2: {
                range: [-xPad, N - 1 + xPad], showticklabels: false, showgrid: false, zeroline: false,
                domain: [xLeft, xRight], anchor: 'y2', showline: false, fixedrange: true
            },
            yaxis2: {
                domain: [hitBot, hitTop], anchor: 'x2',
                range: [0, 1],
                showticklabels: false, showgrid: false, zeroline: false, fixedrange: true
            },
            // Ranked metric panel
            xaxis3: {
                range: [-xPad, N - 1 + xPad],
                domain: [xLeft, xRight], anchor: 'y3',
                showgrid: false,
                tickfont: { size: esMetTickFont.visible !== false ? esMetTickFont.size : 0, family: esMetTickFont.family || fontFam },
                side: 'bottom',
                fixedrange: true
            },
            yaxis3: {
                domain: [metBot, metTop], anchor: 'x3',
                gridcolor: '#eee', gridwidth: 1,
                zeroline: true, zerolinecolor: '#333', zerolinewidth: 0.8,
                tickfont: { size: esMetTickFont.visible !== false ? esMetTickFont.size : 0, family: esMetTickFont.family || fontFam },
                range: [metYMin, metYMax],
                fixedrange: true
            },
            height: s.esPlotHeight || 580,
            margin: { l: 75, r: 15, t: 45, b: 80 },
            font: { family: fontFam },
            paper_bgcolor: s.transparentBg ? 'rgba(0,0,0,0)' : '#fff',
            plot_bgcolor: '#fff',
            shapes: panelShapes,
            annotations: annotations,
            bargap: 0,
            dragmode: false
        };

        // Apply custom dimensions
        const dims = this.getPlotDimensions('esPlot', undefined, layout.height);
        if (dims.width) layout.width = dims.width;
        if (dims.height) layout.height = dims.height;

        // Preserve user-edited annotation positions from previous render
        const plotDiv = document.getElementById('esPlot');
        if (plotDiv && plotDiv.layout && plotDiv.layout.annotations) {
            const prevAnns = plotDiv.layout.annotations;
            // Match annotations by text content and restore x/y positions
            for (const ann of annotations) {
                const match = prevAnns.find(p => p.text === ann.text);
                if (match) {
                    ann.x = match.x;
                    ann.y = match.y;
                    if (match.ax !== undefined) ann.ax = match.ax;
                    if (match.ay !== undefined) ann.ay = match.ay;
                }
            }
        }

        Plotly.newPlot('esPlot',
            [esLine, esZero, rugAnchor, ...traces3, metZero],
            layout,
            {
                responsive: true,
                displaylogo: false,
                modeBarButtonsToRemove: ['lasso2d', 'select2d', 'zoom2d', 'pan2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'resetScale2d'],
                scrollZoom: false,
                doubleClick: false,
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
            sel.dispatchEvent(new Event('change'));
        }
    }

    _initGeneSetSearch() {
        const input = document.getElementById('geneSetSearchInput');
        const sel = document.getElementById('geneSetSelector');
        if (!input || !sel) return;

        // Remove old listeners by replacing element
        const newInput = input.cloneNode(true);
        input.parentNode.replaceChild(newInput, input);

        let debounceTimer = null;
        newInput.addEventListener('input', () => {
            clearTimeout(debounceTimer);
            debounceTimer = setTimeout(() => {
                const query = newInput.value.trim();
                this.populateGeneSetSelector(query || undefined);
            }, 150);
        });

        // Clear search on Escape
        newInput.addEventListener('keydown', (e) => {
            if (e.key === 'Escape') {
                newInput.value = '';
                this.populateGeneSetSelector();
                newInput.blur();
            }
        });
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
        const isUp = result.nes > 0;
        const dirColor = isUp ? '#dc2626' : '#2563eb';
        const metricLabel = document.getElementById('metricColumn').value || 'metric';
        const dataType = this.settings.dataType || 'expression';

        // Build data-type-aware interpretation
        let dirLabel, dirExplanation;
        if (dataType === 'crispr') {
            dirLabel = isUp ? 'Enriched (positive selection)' : 'Depleted (negative selection)';
            dirExplanation = isUp
                ? `Genes in this set tend to have <b>high ${metricLabel}</b> values (enriched in screen). These genes are <b>positively selected</b> — they score high in your screen phenotype (e.g. essential for survival, enriched in a sorted population, etc.).`
                : `Genes in this set tend to have <b>low ${metricLabel}</b> values (depleted in screen). These genes are <b>negatively selected</b> — they score low in your screen phenotype (e.g. dispensable, depleted from a sorted population, etc.).`;
        } else {
            dirLabel = isUp ? 'Upregulated' : 'Downregulated';
            dirExplanation = isUp
                ? `Genes in the gene set tend to have <b>high ${metricLabel}</b> values in your data, meaning this pathway is <b>activated / upregulated</b> in your condition.`
                : `Genes in the gene set tend to have <b>low ${metricLabel}</b> values in your data, meaning this pathway is <b>suppressed / downregulated</b> in your condition.`;
        }

        // FDR interpretation
        let fdrNote = '';
        if (result.fdr < 0.05) fdrNote = '<span style="color:#15803d; font-weight:600;">Significant (FDR < 0.05)</span>';
        else if (result.fdr < 0.25) fdrNote = '<span style="color:#d97706; font-weight:600;">Suggestive (FDR < 0.25)</span>';
        else fdrNote = '<span style="color:#6b7280;">Not significant (FDR ≥ 0.25)</span>';

        const msigdbUrl = `https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/${encodeURIComponent(geneSetName)}.html`;
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
                <span>${result.size} genes in the gene set</span>
                <span style="color: var(--gray-500);">Leading Edge:</span>
                <span>${(result.leadingEdge || []).length} genes</span>
            </div>
            <div style="margin-top: 6px;">
                <a href="${msigdbUrl}" target="_blank" rel="noopener" style="font-size: 0.85em; color: var(--green-600); text-decoration: none; font-weight: 600;">🔗 View on MSigDB</a>
            </div>
            <div style="margin-top: 8px; padding: 6px 8px; background: ${isUp ? '#fef2f2' : '#eff6ff'}; border-radius: 4px; border-left: 3px solid ${dirColor};">
                <div style="font-weight: 600; color: ${dirColor}; margin-bottom: 2px;">${dirLabel}</div>
                <div style="font-size: 0.88em; color: var(--text-secondary); line-height: 1.4;">${dirExplanation}</div>
                <div style="font-size: 0.88em; margin-top: 3px;">${fdrNote}</div>
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
        const collapseThresh = parseInt(document.getElementById('overlapRedundancyFilter')?.value) || 0;
        const fontFam = this.settings.fontFamily + ', sans-serif';

        // Get gene sets filtered by FDR/pval/direction
        const overlapFdrVal = document.getElementById('overlapFdrFilter')?.value || '0.25';
        const overlapPvalVal = document.getElementById('overlapPvalFilter')?.value || '1';
        const overlapDirVal = document.getElementById('overlapDirectionFilter')?.value || 'all';
        const overlapFdrThresh = overlapFdrVal === 'all' ? Infinity : parseFloat(overlapFdrVal);
        const overlapPvalThresh = parseFloat(overlapPvalVal);
        let sigSets = this.results
            .filter(r => {
                if (overlapFdrThresh !== Infinity && r.fdr >= overlapFdrThresh) return false;
                if (overlapPvalThresh < 1 && r.pvalue >= overlapPvalThresh) return false;
                if (overlapDirVal === 'up' && r.nes <= 0) return false;
                if (overlapDirVal === 'down' && r.nes > 0) return false;
                // Respect hidden sets from Select Sets
                if (this._hiddenSets.has(r.name)) return false;
                return true;
            })
            .sort((a, b) => Math.abs(b.nes) - Math.abs(a.nes));

        if (sigSets.length === 0) {
            Plotly.newPlot('overlapHeatmap', [], {
                annotations: [{
                    text: `No gene sets with FDR < ${overlapFdrVal === 'all' ? '∞' : overlapFdrVal} to show`,
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

        // Collapse highly overlapping sets (greedy: keep higher-|NES| representative)
        {
            const kept = [];
            const collapsed = {};
            for (const r of sigSets) {
                if (!setGenes[r.name]) continue;
                let shouldCollapse = false;
                if (collapseThresh > 0) for (const keptSet of kept) {
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
        }

        // Cap at maxSets
        sigSets = sigSets.slice(0, maxSets);

        // Sort alphabetically by default
        sigSets.sort((a, b) => a.name.localeCompare(b.name));

        // Track which sets are visible in the heatmap (for "Send to Results Table")
        this._overlapVisibleSets = new Set(sigSets.map(r => r.name));

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
                    const union = genesI.length + genesJ.length - overlap;
                    const jaccard = union > 0 ? overlap / union : 0;
                    row.push(jaccard);
                    textRow.push(`${overlap} shared<br>J=${jaccard.toFixed(2)}`);
                }
            }
            matrix.push(row);
            textMatrix.push(textRow);
        }

        // Get overlap colorscale from settings
        const overlapColorscales = {
            green: [[0,'#ffffff'],[0.2,'#e0f2e9'],[0.4,'#a8d89a'],[0.6,'#5a9f4a'],[0.8,'#3a7333'],[1.0,'#2d5a27']],
            blues: [[0,'#ffffff'],[0.2,'#deebf7'],[0.4,'#9ecae1'],[0.6,'#4292c6'],[0.8,'#2171b5'],[1.0,'#084594']],
            purples: [[0,'#ffffff'],[0.2,'#efedf5'],[0.4,'#bcbddc'],[0.6,'#807dba'],[0.8,'#6a51a3'],[1.0,'#3f007d']],
            YlGnBu: [[0,'#ffffff'],[0.2,'#edf8b1'],[0.4,'#7fcdbb'],[0.6,'#41b6c4'],[0.8,'#1d91c0'],[1.0,'#0c2c84']]
        };
        const overlapCS = overlapColorscales[this.settings.overlapColorScheme] || overlapColorscales.green;

        const trace = {
            z: matrix,
            x: labels,
            y: labels,
            type: 'heatmap',
            colorscale: overlapCS,
            text: textMatrix,
            hovertemplate: '%{y} vs %{x}<br>%{text}<extra></extra>',
            showscale: true,
            colorbar: {
                title: { text: this._getTextFont('overlap', 'colorbarTitle').visible !== false ? this._getTextFont('overlap', 'colorbarTitle').wrap(this._getTextFont('overlap', 'colorbarTitle').text || 'Jaccard Index') : '', font: { size: this._getTextFont('overlap', 'colorbarTitle').size, family: this._getTextFont('overlap', 'colorbarTitle').family } },
                thickness: 12, len: 0.5
            }
        };

        const overlapTickFont = this._getTextFont('overlap', 'tickFont');
        const overlapTitleFont = this._getTextFont('overlap', 'title');
        const autoSize = Math.max(450, n * 28 + 150);
        const overlapW = parseInt(document.getElementById('overlapPlotWidth')?.value) || 0;
        const overlapH = parseInt(document.getElementById('overlapPlotHeight')?.value) || 0;
        const overlapHeight = overlapH > 0 ? overlapH : autoSize;
        const overlapWidth = overlapW > 0 ? overlapW : autoSize + 50;

        const overlapAnnotations = [];
        if (overlapTitleFont.visible !== false) {
            overlapAnnotations.push({
                text: overlapTitleFont.wrap(overlapTitleFont.text || 'Gene Set Overlap'),
                xref: 'paper', yref: 'paper', x: 0.5, y: 1.02,
                showarrow: false,
                font: { size: overlapTitleFont.size, family: overlapTitleFont.family },
                xanchor: 'center', yanchor: 'bottom'
            });
        }

        const layout = {
            height: overlapHeight,
            width: overlapWidth,
            margin: { l: 10, r: 50, t: 55, b: 10 },
            xaxis: { tickfont: { size: overlapTickFont.visible !== false ? overlapTickFont.size : 0, family: overlapTickFont.family }, tickangle: -45, automargin: true, showgrid: false, fixedrange: true },
            yaxis: { tickfont: { size: overlapTickFont.visible !== false ? overlapTickFont.size : 0, family: overlapTickFont.family }, automargin: true, showgrid: false, autorange: 'reversed', fixedrange: true },
            annotations: overlapAnnotations,
            font: { family: fontFam },
            paper_bgcolor: '#fff',
            plot_bgcolor: '#fff',
            dragmode: false
        };

        Plotly.newPlot('overlapHeatmap', [trace], layout, {
            responsive: true, displaylogo: false,
            modeBarButtonsToRemove: ['lasso2d', 'select2d', 'zoom2d', 'pan2d', 'zoomIn2d', 'zoomOut2d', 'autoScale2d', 'resetScale2d'],
            scrollZoom: false,
            doubleClick: false,
            edits: { annotationPosition: true, annotationText: true, colorbarPosition: true }
        });

        // Hide legacy condensed info panel if it exists
        const condensedEl = document.getElementById('condensedInfo');
        if (condensedEl) condensedEl.style.display = 'none';
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
    // Overlap-based redundancy filtering
    // --------------------------------------------------------

    /** Get the tier (priority) of a gene set: 1=Hallmark, 2=C2/C5, 3=C3/C6, 4=C7/C8, 5=custom */
    _getSetTier(name) {
        if (this.geneSets['hallmark'] && this.geneSets['hallmark'][name]) return 1;
        for (const id of ['c2kegg', 'c2reactome', 'c2wp', 'c2biocarta', 'c2pid', 'c2cgp']) {
            if (this.geneSets[id] && this.geneSets[id][name]) return 2;
        }
        for (const id of ['c5bp', 'c5cc', 'c5mf', 'c5hpo']) {
            if (this.geneSets[id] && this.geneSets[id][name]) return 2;
        }
        if (this.geneSets['c3'] && this.geneSets['c3'][name]) return 3;
        if (this.geneSets['c6'] && this.geneSets['c6'][name]) return 3;
        if (this.geneSets['c7'] && this.geneSets['c7'][name]) return 4;
        if (this.geneSets['c8'] && this.geneSets['c8'][name]) return 4;
        return 5;
    }

    _getSetCollection(name) {
        const colls = [
            ['hallmark', 'H'],
            ['c2kegg', 'KEGG'], ['c2reactome', 'Reactome'], ['c2wp', 'WikiPathways'],
            ['c2biocarta', 'BioCarta'], ['c2pid', 'PID'], ['c2cgp', 'CGP'],
            ['c5bp', 'GO:BP'], ['c5cc', 'GO:CC'], ['c5mf', 'GO:MF'], ['c5hpo', 'HPO'],
            ['c3', 'C3'], ['c6', 'C6'], ['c7', 'C7'], ['c8', 'C8']
        ];
        for (const [id, label] of colls) {
            if (this.geneSets[id] && this.geneSets[id][name]) return label;
        }
        return 'Custom';
    }

    /** Build pairwise overlap cache for all results. Called after GSEA completes. */
    _buildOverlapCache() {
        if (!this.results || !this.rankedList) return;
        const activeGeneSets = this.getActiveGeneSets();
        if (!activeGeneSets) return;

        const rankedGenesUpper = new Set(this.rankedList.genes.map(g => g.toUpperCase()));
        const setGenes = {};

        // For large result sets, only compute overlap for top results by FDR
        const MAX_OVERLAP_SETS = 200;
        let overlapResults = this.results;
        this._overlapCapped = false;
        if (overlapResults.length > MAX_OVERLAP_SETS) {
            overlapResults = [...overlapResults]
                .sort((a, b) => a.fdr - b.fdr)
                .slice(0, MAX_OVERLAP_SETS);
            this._overlapCapped = true;
        }

        for (const r of overlapResults) {
            if (activeGeneSets[r.name]) {
                setGenes[r.name] = activeGeneSets[r.name]
                    .map(g => g.toUpperCase())
                    .filter(g => rankedGenesUpper.has(g));
            }
        }

        // Build pairwise Jaccard map: key = "nameA|||nameB" (sorted), value = jaccard
        const pairwise = new Map();
        const names = overlapResults.map(r => r.name).filter(n => setGenes[n]);
        for (let i = 0; i < names.length; i++) {
            for (let j = i + 1; j < names.length; j++) {
                const a = names[i], b = names[j];
                const genesA = setGenes[a], genesB = setGenes[b];
                if (!genesA || !genesB) continue;
                const inter = this._computeOverlap(genesA, genesB);
                const union = genesA.length + genesB.length - inter;
                const jaccard = union > 0 ? inter / union : 0;
                if (jaccard > 0.05) { // only store non-trivial overlaps
                    const key = a < b ? `${a}|||${b}` : `${b}|||${a}`;
                    pairwise.set(key, jaccard);
                }
            }
        }

        this._overlapCache = { pairwise, setGenes };
    }

    /** Get Jaccard overlap between two gene sets from cache */
    _getCachedJaccard(nameA, nameB) {
        if (!this._overlapCache) return 0;
        const key = nameA < nameB ? `${nameA}|||${nameB}` : `${nameB}|||${nameA}`;
        return this._overlapCache.pairwise.get(key) || 0;
    }

    /**
     * Apply overlap filter: given a list of results, remove redundant sets.
     * For each pair with Jaccard > threshold, keep the one with higher tier (or higher |NES| if same tier).
     * Returns filtered array.
     */
    _applyOverlapFilter(results, thresholdPct) {
        if (!thresholdPct || thresholdPct <= 0 || !this._overlapCache) return results;
        const threshold = thresholdPct / 100;

        // Sort by tier (ascending = better first), then |NES| descending
        const sorted = results.slice().sort((a, b) => {
            const tA = this._getSetTier(a.name), tB = this._getSetTier(b.name);
            if (tA !== tB) return tA - tB;
            return Math.abs(b.nes) - Math.abs(a.nes);
        });

        const kept = [];
        const removed = new Set();

        for (const r of sorted) {
            if (removed.has(r.name)) continue;
            if (this._hiddenSets.has(r.name)) continue;
            kept.push(r);
            // Mark all lower-priority sets overlapping with this one
            for (const other of sorted) {
                if (other.name === r.name || removed.has(other.name)) continue;
                if (kept.some(k => k.name === other.name)) continue;
                const jaccard = this._getCachedJaccard(r.name, other.name);
                if (jaccard >= threshold) {
                    removed.add(other.name);
                }
            }
        }

        return kept;
    }

    /** Compute overlap clusters: connected components at a given Jaccard threshold */
    _computeOverlapClusters(results, threshold) {
        if (!this._overlapCache || threshold <= 0) return {};
        const clusters = {}; // name → cluster id
        let nextCluster = 1;
        const sorted = results.slice().sort((a, b) => {
            const tA = this._getSetTier(a.name), tB = this._getSetTier(b.name);
            if (tA !== tB) return tA - tB;
            return Math.abs(b.nes) - Math.abs(a.nes);
        });
        for (const r of sorted) {
            if (clusters[r.name]) continue;
            // Find all sets overlapping with this one (transitively via BFS)
            const queue = [r.name];
            const visited = new Set([r.name]);
            while (queue.length > 0) {
                const current = queue.shift();
                for (const other of sorted) {
                    if (visited.has(other.name)) continue;
                    const j = this._getCachedJaccard(current, other.name);
                    if (j >= threshold) {
                        visited.add(other.name);
                        queue.push(other.name);
                    }
                }
            }
            if (visited.size > 1) {
                for (const name of visited) {
                    clusters[name] = nextCluster;
                }
                nextCluster++;
            }
        }
        return clusters;
    }

    /** Get the best representative per overlap cluster */
    _getClusterRepresentatives(clusters) {
        const reps = new Set();
        const clusterSets = {};
        for (const [name, cid] of Object.entries(clusters)) {
            if (!clusterSets[cid]) clusterSets[cid] = [];
            clusterSets[cid].push(name);
        }
        for (const members of Object.values(clusterSets)) {
            // Sort by tier (ascending) then |NES| (descending) — best is first
            members.sort((a, b) => {
                const tA = this._getSetTier(a), tB = this._getSetTier(b);
                if (tA !== tB) return tA - tB;
                const rA = this.results.find(r => r.name === a);
                const rB = this.results.find(r => r.name === b);
                return Math.abs((rB?.nes || 0)) - Math.abs((rA?.nes || 0));
            });
            reps.add(members[0]);
        }
        return reps;
    }

    /** Open the gene set selector popup for choosing which sets to show/hide */
    openGeneSetFilter(context) {
        // context: 'bubble' or 'table'
        if (!this.results) return;
        this._gsfContext = context;
        this._gsfSortCol = this._gsfSortCol || 'nes';
        this._gsfSortAsc = this._gsfSortAsc !== undefined ? this._gsfSortAsc : false;

        // Set defaults on first open: FDR < 0.25, overlap > 30%, auto-select
        if (!this._gsfInitialized) {
            this._gsfFdrFilter = '0.25';
            this._gsfClusterThreshold = 0.3;
            this._gsfInitialized = true;
        }

        this._renderGeneSetFilter();

        // Auto-select best per cluster on first open
        if (this._gsfAutoSelectOnOpen) return; // already done
        this._gsfAutoSelectOnOpen = true;

        // First: apply FDR filter to _hiddenSets — hide sets that don't pass
        const fdrT = this._gsfFdrFilter && this._gsfFdrFilter !== 'all'
            ? parseFloat(this._gsfFdrFilter) : Infinity;
        for (const r of this.results) {
            if (fdrT < Infinity && r.fdr >= fdrT) {
                this._hiddenSets.add(r.name);
            } else {
                this._hiddenSets.delete(r.name);
            }
        }

        // Then: among passing sets, auto-select best per overlap cluster
        const thresh = this._gsfClusterThreshold || 0.3;
        if (thresh > 0) {
            const clusterInput = this.results.filter(r => r.fdr < fdrT);
            const cappedInput = clusterInput.length > 500
                ? clusterInput.slice().sort((a, b) => Math.abs(b.nes) - Math.abs(a.nes)).slice(0, 500)
                : clusterInput;
            const cls = this._computeOverlapClusters(cappedInput, thresh);
            const reps = this._getClusterRepresentatives(cls);
            for (const [name, cid] of Object.entries(cls)) {
                if (reps.has(name)) {
                    this._hiddenSets.delete(name);
                } else {
                    this._hiddenSets.add(name);
                }
            }
            this._renderGeneSetFilter();
        }
    }

    _renderGeneSetFilter() {
        const popup = document.getElementById('geneSetFilterPopup');
        const body = document.getElementById('geneSetFilterBody');
        const gsfFdrVal = this._gsfFdrFilter || 'all';
        const gsfPvalVal = this._gsfPvalFilter || 'all';
        const gsfClusterThresh = this._gsfClusterThreshold !== undefined ? this._gsfClusterThreshold : 0.3;

        // Compute clusters (on filtered results to avoid O(n²) on thousands of sets)
        let clusterInput = this.results;
        if (gsfFdrVal !== 'all') {
            const t = parseFloat(gsfFdrVal);
            clusterInput = clusterInput.filter(r => r.fdr < t);
        }
        if (gsfPvalVal !== 'all') {
            const t = parseFloat(gsfPvalVal);
            clusterInput = clusterInput.filter(r => r.pvalue < t);
        }
        if (clusterInput.length > 500) {
            clusterInput = clusterInput.slice().sort((a, b) => Math.abs(b.nes) - Math.abs(a.nes)).slice(0, 500);
        }
        const rawClusters = this._computeOverlapClusters(clusterInput, gsfClusterThresh);
        // Renumber clusters sequentially so displayed sets start from 1
        const clusterRemap = {};
        let nextId = 1;
        // Assign new IDs in order of first appearance (sorted by |NES|)
        const clusterOrder = clusterInput.slice().sort((a, b) => Math.abs(b.nes) - Math.abs(a.nes));
        for (const r of clusterOrder) {
            const oldId = rawClusters[r.name];
            if (oldId && !clusterRemap[oldId]) {
                clusterRemap[oldId] = nextId++;
            }
        }
        const clusters = {};
        for (const [name, oldId] of Object.entries(rawClusters)) {
            clusters[name] = clusterRemap[oldId] || oldId;
        }

        // Sort results
        let sorted = this.results.slice();
        const sortCol = this._gsfSortCol || 'nes';
        const sortAsc = this._gsfSortAsc || false;
        sorted.sort((a, b) => {
            let vA, vB;
            if (sortCol === 'name') { vA = a.name.toLowerCase(); vB = b.name.toLowerCase(); }
            else if (sortCol === 'size') { vA = a.size; vB = b.size; }
            else if (sortCol === 'nes') { vA = Math.abs(a.nes); vB = Math.abs(b.nes); }
            else if (sortCol === 'fdr') { vA = a.fdr; vB = b.fdr; }
            else if (sortCol === 'pval') { vA = a.pvalue; vB = b.pvalue; }
            else if (sortCol === 'cluster') { vA = clusters[a.name] || 999; vB = clusters[b.name] || 999; }
            else { vA = a[sortCol]; vB = b[sortCol]; }
            if (vA < vB) return sortAsc ? -1 : 1;
            if (vA > vB) return sortAsc ? 1 : -1;
            return 0;
        });

        // Filter by FDR
        if (gsfFdrVal !== 'all') {
            const thresh = parseFloat(gsfFdrVal);
            sorted = sorted.filter(r => r.fdr < thresh);
        }
        // Filter by p-value
        if (gsfPvalVal !== 'all') {
            const thresh = parseFloat(gsfPvalVal);
            sorted = sorted.filter(r => r.pvalue < thresh);
        }
        // Filter by search query
        const searchQ = (this._gsfSearchVal || '').toLowerCase();
        if (searchQ) {
            sorted = sorted.filter(r => r.name.toLowerCase().includes(searchQ));
        }

        const clusterColors = ['#d97706', '#7c3aed', '#059669', '#dc2626', '#2563eb', '#db2777', '#ca8a04', '#0891b2'];

        const sortArrow = (col) => {
            if (sortCol === col) return sortAsc ? ' &#9650;' : ' &#9660;';
            return ' <span style="opacity:0.3">&#9650;</span>';
        };
        const thStyle = 'padding: 4px 6px; cursor: pointer; user-select: none; white-space: nowrap;';

        let html = `<div style="display: flex; gap: 8px; margin-bottom: 8px; align-items: center; flex-wrap: wrap;">`;
        html += `<input type="text" id="gsfSearch" class="form-control" placeholder="Filter..." style="flex: 1; max-width: 180px; font-size: 0.85em;" value="${this._gsfSearchVal || ''}">`;
        html += `<select class="form-control" id="gsfFdrFilter" style="width: auto; font-size: 0.85em;">`;
        html += `<option value="all"${gsfFdrVal === 'all' ? ' selected' : ''}>All FDR</option>`;
        html += `<option value="0.5"${gsfFdrVal === '0.5' ? ' selected' : ''}>FDR < 0.5</option>`;
        html += `<option value="0.25"${gsfFdrVal === '0.25' ? ' selected' : ''}>FDR < 0.25</option>`;
        html += `<option value="0.1"${gsfFdrVal === '0.1' ? ' selected' : ''}>FDR < 0.1</option>`;
        html += `<option value="0.05"${gsfFdrVal === '0.05' ? ' selected' : ''}>FDR < 0.05</option>`;
        html += `<option value="0.01"${gsfFdrVal === '0.01' ? ' selected' : ''}>FDR < 0.01</option>`;
        html += `</select>`;
        html += `<select class="form-control" id="gsfPvalFilter" style="width: auto; font-size: 0.85em;">`;
        html += `<option value="all"${gsfPvalVal === 'all' ? ' selected' : ''}>All p-val</option>`;
        html += `<option value="0.1"${gsfPvalVal === '0.1' ? ' selected' : ''}>p < 0.1</option>`;
        html += `<option value="0.05"${gsfPvalVal === '0.05' ? ' selected' : ''}>p < 0.05</option>`;
        html += `<option value="0.01"${gsfPvalVal === '0.01' ? ' selected' : ''}>p < 0.01</option>`;
        html += `<option value="0.001"${gsfPvalVal === '0.001' ? ' selected' : ''}>p < 0.001</option>`;
        html += `</select>`;
        html += `<select class="form-control" id="gsfClusterThresh" style="width: auto; font-size: 0.85em;" title="Jaccard threshold for overlap clusters">`;
        for (const v of [0, 0.1, 0.2, 0.3, 0.5]) {
            html += `<option value="${v}"${gsfClusterThresh === v ? ' selected' : ''}>${v === 0 ? 'No clusters' : `Overlap > ${(v * 100).toFixed(0)}%`}</option>`;
        }
        html += `</select>`;
        html += `<button class="btn btn-outline btn-sm" id="gsfAutoSelect" title="Auto-select best representative per overlap cluster">Auto-select</button>`;
        html += `<button class="btn btn-outline btn-sm" id="gsfShowAll">Show all</button>`;
        html += `<button class="btn btn-outline btn-sm" id="gsfHideAll">Hide all</button>`;
        html += `<button class="btn btn-sm" id="gsfViewInOverview" style="background: var(--green-50); border: 1px solid var(--green-600); color: var(--green-700);" title="Pin visible sets and view in Overview lollipop plot">📊 View in Overview</button>`;
        html += `</div>`;

        html += `<div style="max-height: 400px; overflow-y: auto; border: 1px solid #eee; border-radius: 4px;">`;
        html += `<table style="width: 100%; border-collapse: collapse; font-size: 0.82em;">`;
        html += `<thead><tr style="background: #f9fafb; position: sticky; top: 0; z-index: 1;">`;
        html += `<th style="padding: 4px 6px; text-align: left; width: 30px;"></th>`;
        html += `<th class="gsf-sort" data-gsf-col="name" style="${thStyle} text-align: left;">Gene Set${sortArrow('name')}</th>`;
        html += `<th class="gsf-sort" data-gsf-col="size" style="${thStyle} text-align: right; width: 45px;">Size${sortArrow('size')}</th>`;
        html += `<th class="gsf-sort" data-gsf-col="nes" style="${thStyle} text-align: right; width: 50px;">NES${sortArrow('nes')}</th>`;
        html += `<th class="gsf-sort" data-gsf-col="fdr" style="${thStyle} text-align: right; width: 55px;">FDR${sortArrow('fdr')}</th>`;
        html += `<th class="gsf-sort" data-gsf-col="pval" style="${thStyle} text-align: right; width: 55px;">p-val${sortArrow('pval')}</th>`;
        html += `<th style="padding: 4px 6px; text-align: center; width: 35px;">Coll.</th>`;
        if (gsfClusterThresh > 0) html += `<th class="gsf-sort" data-gsf-col="cluster" style="${thStyle} text-align: center; width: 55px;">Cluster${sortArrow('cluster')}</th>`;
        html += `</tr></thead><tbody>`;

        const MAX_DISPLAY_ROWS = 500;
        let displayCount = 0;
        const totalMatching = sorted.length;
        for (const r of sorted) {
            if (displayCount >= MAX_DISPLAY_ROWS) break;
            const hidden = this._hiddenSets.has(r.name);
            const coll = this._getSetCollection(r.name);
            const nesColor = r.nes > 0 ? '#dc2626' : '#2563eb';
            const fdrStr = r.fdr < 0.001 ? r.fdr.toExponential(1) : r.fdr.toFixed(3);
            const cid = clusters[r.name];
            const clusterBadge = (cid && gsfClusterThresh > 0)
                ? `<span style="display:inline-block; min-width:18px; text-align:center; padding:1px 5px; border-radius:10px; font-size:0.8em; font-weight:600; color:white; background:${clusterColors[(cid - 1) % clusterColors.length]}">${cid}</span>`
                : '';

            html += `<tr class="gsf-row" data-name="${this._escapeAttr(r.name)}" style="border-bottom: 1px solid #f5f5f5; ${hidden ? 'opacity: 0.45;' : ''}">`;
            html += `<td style="padding: 3px 6px;"><input type="checkbox" class="gsf-check" data-name="${this._escapeAttr(r.name)}" ${hidden ? '' : 'checked'}></td>`;
            html += `<td style="padding: 3px 6px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; max-width: 220px;" title="${this._escapeAttr(r.name)}">${this.cleanName(r.name)}</td>`;
            html += `<td style="padding: 3px 6px; text-align: right; font-family: Roboto Mono, monospace; color: var(--gray-500);">${r.size}</td>`;
            html += `<td style="padding: 3px 6px; text-align: right; color: ${nesColor}; font-family: Roboto Mono, monospace;">${r.nes.toFixed(2)}</td>`;
            html += `<td style="padding: 3px 6px; text-align: right; font-family: Roboto Mono, monospace;">${fdrStr}</td>`;
            const pvalStr = r.pvalue < 0.001 ? r.pvalue.toExponential(1) : r.pvalue.toFixed(3);
            html += `<td style="padding: 3px 6px; text-align: right; font-family: Roboto Mono, monospace;">${pvalStr}</td>`;
            html += `<td style="padding: 3px 6px; text-align: center; font-size: 0.85em; color: var(--gray-500);">${coll}</td>`;
            if (gsfClusterThresh > 0) html += `<td style="padding: 3px 6px; text-align: center;">${clusterBadge}</td>`;
            html += `</tr>`;
            displayCount++;
        }

        if (totalMatching > MAX_DISPLAY_ROWS) {
            html += `<tr><td colspan="8" style="padding: 8px; text-align: center; color: var(--gray-500); font-size: 0.85em;">Showing ${MAX_DISPLAY_ROWS} of ${totalMatching} matching sets. Use the search box to narrow results.</td></tr>`;
        }
        html += `</tbody></table></div>`;
        html += `<div style="margin-top: 8px; display: flex; gap: 6px; justify-content: flex-end;">`;
        const visibleCount = this.results.filter(r => !this._hiddenSets.has(r.name)).length;
        html += `<span style="flex: 1; font-size: 0.82em; color: var(--gray-500); align-self: center;" id="gsfCount">${visibleCount} of ${this.results.length} selected</span>`;
        html += `<button class="btn btn-outline btn-sm" id="gsfCancel">Cancel</button>`;
        html += `<button class="btn btn-sm" id="gsfApply" style="background: var(--green-600); color: white;">Apply</button>`;
        html += `</div>`;

        body.innerHTML = html;
        popup.style.display = 'block';
        document.getElementById('geneSetFilterBackdrop').style.display = 'block';

        // Restore focus to search input if user was typing
        if (searchQ) {
            const searchInput = document.getElementById('gsfSearch');
            if (searchInput) {
                searchInput.focus();
                searchInput.setSelectionRange(searchInput.value.length, searchInput.value.length);
            }
        }

        // Wire up events (re-bind each time since we rebuild the entire DOM)
        // Use event delegation on body
        const self = this;
        body.onchange = function(e) {
            if (e.target.classList.contains('gsf-check')) {
                const name = e.target.dataset.name;
                const row = e.target.closest('.gsf-row');
                if (e.target.checked) {
                    self._hiddenSets.delete(name);
                    if (row) row.style.opacity = '1';
                } else {
                    self._hiddenSets.add(name);
                    if (row) row.style.opacity = '0.45';
                }
                const total = self.results.length;
                document.getElementById('gsfCount').textContent = `${total - self._hiddenSets.size} of ${total} selected`;
            }
            if (e.target.id === 'gsfFdrFilter') {
                self._gsfFdrFilter = e.target.value;
                // Auto-update selection: show sets passing filter, hide those that don't
                const fdrT = self._gsfFdrFilter === 'all' ? Infinity : parseFloat(self._gsfFdrFilter);
                for (const r of self.results) {
                    if (r.fdr < fdrT) {
                        self._hiddenSets.delete(r.name);
                    } else {
                        self._hiddenSets.add(r.name);
                    }
                }
                // Also apply p-value filter if set
                const pT = self._gsfPvalFilter && self._gsfPvalFilter !== 'all' ? parseFloat(self._gsfPvalFilter) : Infinity;
                if (pT < Infinity) {
                    for (const r of self.results) {
                        if (r.pvalue >= pT) self._hiddenSets.add(r.name);
                    }
                }
                self._renderGeneSetFilter();
            }
            if (e.target.id === 'gsfPvalFilter') {
                self._gsfPvalFilter = e.target.value;
                // Auto-update selection: show sets passing filter, hide those that don't
                const pT = self._gsfPvalFilter === 'all' ? Infinity : parseFloat(self._gsfPvalFilter);
                for (const r of self.results) {
                    if (r.pvalue < pT) {
                        self._hiddenSets.delete(r.name);
                    } else {
                        self._hiddenSets.add(r.name);
                    }
                }
                // Also apply FDR filter if set
                const fdrT = self._gsfFdrFilter && self._gsfFdrFilter !== 'all' ? parseFloat(self._gsfFdrFilter) : Infinity;
                if (fdrT < Infinity) {
                    for (const r of self.results) {
                        if (r.fdr >= fdrT) self._hiddenSets.add(r.name);
                    }
                }
                self._renderGeneSetFilter();
            }
            if (e.target.id === 'gsfClusterThresh') {
                self._gsfClusterThreshold = parseFloat(e.target.value);
                self._renderGeneSetFilter();
            }
        };
        body.onclick = function(e) {
            const sortTh = e.target.closest('.gsf-sort');
            if (sortTh) {
                const col = sortTh.dataset.gsfCol;
                if (self._gsfSortCol === col) {
                    self._gsfSortAsc = !self._gsfSortAsc;
                } else {
                    self._gsfSortCol = col;
                    self._gsfSortAsc = col === 'name' || col === 'fdr' || col === 'pval';
                }
                self._renderGeneSetFilter();
                return;
            }
            if (e.target.id === 'gsfApply') {
                popup.style.display = 'none';
                document.getElementById('geneSetFilterBackdrop').style.display = 'none';
                // Sync Select Sets filters to all tabs
                const syncFdr = self._gsfFdrFilter || 'all';
                const syncPval = self._gsfPvalFilter || 'all';
                // Map Select Sets FDR/pval values to per-tab dropdown values
                const setVal = (id, val) => { const el = document.getElementById(id); if (el) el.value = val; };
                // Results Table (uses 'all' for no filter)
                setVal('fdrFilter', syncFdr === 'all' ? 'all' : syncFdr);
                setVal('pvalueFilter', syncPval === 'all' ? 'all' : syncPval);
                // Overview (uses '1' for no filter)
                setVal('overviewFdrFilter', syncFdr === 'all' ? '1' : syncFdr);
                setVal('overviewPvalFilter', syncPval === 'all' ? '1' : syncPval);
                // Enrichment Score Plot (uses '1' for no filter)
                setVal('esFdrFilter', syncFdr === 'all' ? '1' : syncFdr);
                setVal('esPvalFilter', syncPval === 'all' ? '1' : syncPval);
                // If Select Sets used overlap clustering, reset per-tab overlap filters
                if (self._gsfClusterThreshold > 0) {
                    setVal('tableOverlapFilter', '0');
                    setVal('overviewOverlapFilter', '0');
                    self._overlapFilterThreshold = 0;
                }
                // Apply to ALL tabs — hidden sets are shared
                self.renderBubblePlot();
                self.filterAndRenderTable();
                if (self.results) self.renderOverlapHeatmap();
                // Also refresh enrichment plot gene set selector
                self.populateGeneSetSelector();
                self._initGeneSetSearch();
            } else if (e.target.id === 'gsfCancel') {
                popup.style.display = 'none';
                document.getElementById('geneSetFilterBackdrop').style.display = 'none';
            } else if (e.target.id === 'gsfShowAll') {
                self._hiddenSets.clear();
                body.querySelectorAll('.gsf-check').forEach(c => { c.checked = true; });
                body.querySelectorAll('.gsf-row').forEach(r => { r.style.opacity = '1'; });
                document.getElementById('gsfCount').textContent = `${self.results.length} of ${self.results.length} selected`;
            } else if (e.target.id === 'gsfHideAll') {
                self.results.forEach(r => self._hiddenSets.add(r.name));
                body.querySelectorAll('.gsf-check').forEach(c => { c.checked = false; });
                body.querySelectorAll('.gsf-row').forEach(r => { r.style.opacity = '0.45'; });
                document.getElementById('gsfCount').textContent = `0 of ${self.results.length} selected`;
            } else if (e.target.id === 'gsfAutoSelect') {
                // Auto-select: show best per cluster, hide others in clusters
                const thresh = self._gsfClusterThreshold !== undefined ? self._gsfClusterThreshold : 0.3;
                if (thresh <= 0) return;
                const cls = self._computeOverlapClusters(self.results, thresh);
                const reps = self._getClusterRepresentatives(cls);
                // For each gene set in a cluster: show rep, hide non-rep
                for (const [name, cid] of Object.entries(cls)) {
                    if (reps.has(name)) {
                        self._hiddenSets.delete(name);
                    } else {
                        self._hiddenSets.add(name);
                    }
                }
                self._renderGeneSetFilter();
            } else if (e.target.id === 'gsfViewInOverview') {
                // Pin all currently checked sets and show ONLY them in Overview
                self._pinnedBubbleSets.clear();
                body.querySelectorAll('.gsf-check:checked').forEach(cb => {
                    self._pinnedBubbleSets.add(cb.dataset.name);
                });
                popup.style.display = 'none';
                document.getElementById('geneSetFilterBackdrop').style.display = 'none';
                self._showOnlyPinned = true;
                self.renderBubblePlot();
                self.filterAndRenderTable();
                self.showTab('overview');
            }
        };
        let gsfSearchTimer = null;
        body.oninput = function(e) {
            if (e.target.id === 'gsfSearch') {
                self._gsfSearchVal = e.target.value;
                clearTimeout(gsfSearchTimer);
                gsfSearchTimer = setTimeout(() => self._renderGeneSetFilter(), 200);
            }
        };
    }

    // --------------------------------------------------------
    // Results Table
    // --------------------------------------------------------
    filterAndRenderTable() {
        if (!this.results) return;

        const query = document.getElementById('tableFilter').value.toLowerCase();
        const fdrVal = document.getElementById('fdrFilter').value;
        const pvalVal = document.getElementById('pvalueFilter').value;
        const dirVal = document.getElementById('directionFilter').value;

        let filtered = this.results.slice();

        // Apply hidden sets filter
        if (this._hiddenSets.size > 0) {
            filtered = filtered.filter(r => !this._hiddenSets.has(r.name));
        }
        // Apply overlap redundancy filter
        const overlapVal = document.getElementById('tableOverlapFilter') ? document.getElementById('tableOverlapFilter').value : '0';
        if (overlapVal !== '0') {
            filtered = this._applyOverlapFilter(filtered, parseInt(overlapVal));
        }

        if (query) {
            filtered = filtered.filter(r => r.name.toLowerCase().includes(query));
        }
        if (fdrVal !== 'all') {
            const thresh = parseFloat(fdrVal);
            filtered = filtered.filter(r => r.fdr < thresh);
        }
        if (pvalVal !== 'all') {
            const thresh = parseFloat(pvalVal);
            filtered = filtered.filter(r => r.pvalue < thresh);
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
            const isPinned = this._pinnedBubbleSets.has(r.name);
            tr.innerHTML = `
                <td style="width:28px;padding:2px;text-align:center;cursor:pointer" class="pin-cell" title="${isPinned ? 'Unpin from overview' : 'Pin to overview (lollipop plot)'}">
                    <span style="opacity:${isPinned ? '1' : '0.25'};font-size:13px">&#128204;</span>
                </td>
                <td title="${r.name}">${this.cleanName(r.name)}</td>
                <td>${r.size}</td>
                <td class="${r.es >= 0 ? 'positive' : 'negative'}">${r.es.toFixed(4)}</td>
                <td class="${r.nes >= 0 ? 'positive' : 'negative'}"><strong>${r.nes.toFixed(3)}</strong></td>
                <td>${this.formatPval(r.pvalue)}</td>
                <td>${this.formatPval(r.fdr)}</td>
                <td class="leading-edge-cell" title="${(r.leadingEdge || []).join(', ')}">${(r.leadingEdge || []).length} genes</td>
            `;
            // Pin/unpin click on first cell
            tr.querySelector('.pin-cell').addEventListener('click', (e) => {
                e.stopPropagation();
                if (this._pinnedBubbleSets.has(r.name)) {
                    this._pinnedBubbleSets.delete(r.name);
                } else {
                    this._pinnedBubbleSets.add(r.name);
                }
                // Update pin icon immediately
                const pinSpan = e.currentTarget.querySelector('span');
                if (pinSpan) pinSpan.style.opacity = this._pinnedBubbleSets.has(r.name) ? '1' : '0.25';
                e.currentTarget.title = this._pinnedBubbleSets.has(r.name) ? 'Unpin from overview' : 'Pin to overview (lollipop plot)';
                // Show/hide "View pinned" button
                const vpBtn = document.getElementById('viewPinnedBtn');
                if (vpBtn) vpBtn.style.display = this._pinnedBubbleSets.size > 0 ? '' : 'none';
                this.renderBubblePlot();
            });
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
            if (!arrow) return; // skip pin column (no sort arrow)
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
        if (document.getElementById('checkC2kegg').checked) collections.push('KEGG (C2)');
        if (document.getElementById('checkC2reactome').checked) collections.push('Reactome (C2)');
        if (document.getElementById('checkC2wp').checked) collections.push('WikiPathways (C2)');
        if (document.getElementById('checkC2biocarta').checked) collections.push('BioCarta (C2)');
        if (document.getElementById('checkC2pid').checked) collections.push('PID (C2)');
        if (document.getElementById('checkC2cgp').checked) collections.push('CGP (C2)');
        if (document.getElementById('checkC5bp').checked) collections.push('GO:BP (C5)');
        if (document.getElementById('checkC5cc').checked) collections.push('GO:CC (C5)');
        if (document.getElementById('checkC5mf').checked) collections.push('GO:MF (C5)');
        if (document.getElementById('checkC5hpo').checked) collections.push('HPO (C5)');
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

        // Read the ACTUAL rendered dimensions from Plotly's internal layout.
        // _fullLayout includes computed automargin adjustments.
        const fullLayout = plotEl._fullLayout || {};
        let exportWidth = fullLayout.width || plotEl.offsetWidth || 800;
        let exportHeight = fullLayout.height || plotEl.offsetHeight || 500;

        // Use Plotly's computed _size (includes automargin) to ensure nothing is cropped.
        // _size.l/.r/.t/.b are the ACTUAL margins after automargin expansion.
        const size = fullLayout._size;
        if (size) {
            const fullW = Math.ceil(size.l + size.w + size.r) + 10;
            const fullH = Math.ceil(size.t + size.h + size.b) + 10;
            exportWidth = Math.max(exportWidth, fullW);
            exportHeight = Math.max(exportHeight, fullH);
        }

        const opts = {
            format: format,
            width: exportWidth,
            height: exportHeight
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

    // --------------------------------------------------------
    // CSV Export (per-plot data)
    // --------------------------------------------------------
    exportPlotCSV(plotType) {
        if (!this.results || !this.rankedList) return;
        this.readSettings();
        let csv = '', filename = '';

        switch (plotType) {
            case 'bubble': {
                const fdrThresh = parseFloat(this.settings.fdrDisplayThreshold);
                const topN = this.settings.topN;
                let filtered = this.results.slice();
                if (fdrThresh < 1) filtered = filtered.filter(r => r.fdr < fdrThresh);
                const top = filtered.sort((a, b) => Math.abs(b.nes) - Math.abs(a.nes)).slice(0, topN);
                csv = 'Gene Set,NES,ES,p-value,FDR,Size,Leading Edge\n';
                csv += top.map(r => [
                    `"${r.name}"`, r.nes.toFixed(4), r.es.toFixed(4),
                    r.pvalue.toExponential(4), r.fdr.toExponential(4), r.size,
                    `"${(r.leadingEdge || []).join(', ')}"`
                ].join(',')).join('\n');
                filename = 'gsea_overview.csv';
                break;
            }
            case 'es': {
                const sel = document.getElementById('geneSetSelector').value;
                const result = this.results.find(r => r.name === sel);
                if (!result) return;
                csv = 'Rank,Gene,Metric,Running_ES,In_Set,Leading_Edge\n';
                const genes = this.rankedList.genes;
                const metrics = this.rankedList.metrics;
                const hitSet = new Set(result.hits);
                const peakIdx = result.runningES.indexOf(result.es);
                for (let i = 0; i < genes.length; i++) {
                    const inSet = hitSet.has(i);
                    const isLE = inSet && (result.es >= 0 ? i <= peakIdx : i >= peakIdx);
                    csv += [i + 1, genes[i], metrics[i].toFixed(6),
                        result.runningES[i].toFixed(6), inSet ? 'Yes' : 'No',
                        isLE ? 'Yes' : 'No'].join(',') + '\n';
                }
                filename = `gsea_es_${sel.replace(/[^a-zA-Z0-9]/g, '_').substring(0, 40)}.csv`;
                break;
            }
            case 'ranked': {
                csv = 'Rank,Gene,Metric\n';
                const genes = this.rankedList.genes;
                const metrics = this.rankedList.metrics;
                for (let i = 0; i < genes.length; i++) {
                    csv += [i + 1, genes[i], metrics[i].toFixed(6)].join(',') + '\n';
                }
                filename = 'gsea_ranked_list.csv';
                break;
            }
            case 'overlap': {
                const maxSets = parseInt(document.getElementById('overlapMaxSets').value) || 20;
                let sigSets = this.results.filter(r => r.fdr < 0.25)
                    .sort((a, b) => Math.abs(b.nes) - Math.abs(a.nes)).slice(0, maxSets);
                const activeGeneSets = this.getActiveGeneSets();
                const setGenes = {};
                const rankedGenesUpper = new Set(this.rankedList.genes.map(g => g.toUpperCase()));
                for (const r of sigSets) {
                    if (activeGeneSets[r.name]) {
                        setGenes[r.name] = activeGeneSets[r.name].map(g => g.toUpperCase()).filter(g => rankedGenesUpper.has(g));
                    }
                }
                sigSets = sigSets.filter(r => setGenes[r.name] && setGenes[r.name].length > 0);
                const names = sigSets.map(r => this.cleanName(r.name));
                csv = ',' + names.join(',') + '\n';
                for (let i = 0; i < sigSets.length; i++) {
                    const row = [names[i]];
                    for (let j = 0; j < sigSets.length; j++) {
                        if (i === j) { row.push('1.0000'); continue; }
                        const setA = new Set(setGenes[sigSets[i].name]);
                        const setB = new Set(setGenes[sigSets[j].name]);
                        let inter = 0;
                        for (const g of setA) if (setB.has(g)) inter++;
                        const union = setA.size + setB.size - inter;
                        row.push(union > 0 ? (inter / union).toFixed(4) : '0.0000');
                    }
                    csv += row.join(',') + '\n';
                }
                filename = 'gsea_overlap_matrix.csv';
                break;
            }
            default: return;
        }

        this.downloadBlob(new Blob([csv], { type: 'text/csv' }), filename);
    }

    // --------------------------------------------------------
    // R Script Export
    // --------------------------------------------------------
    exportRScript(plotType) {
        if (!this.results || !this.rankedList) return;
        this.readSettings();

        let script = '';
        const date = new Date().toISOString().split('T')[0];
        const header = `# ============================================================\n# Generated by Enrich — GSEA Web Tool\n# ${date}\n# https://fredrikwermeling.github.io/enrich/\n# ============================================================\n\n`;

        switch (plotType) {
            case 'bubble': script = header + this._rScriptBubble(); break;
            case 'es':     script = header + this._rScriptES(); break;
            case 'ranked': script = header + this._rScriptRanked(); break;
            case 'overlap': script = header + this._rScriptOverlap(); break;
            default: return;
        }

        const blob = new Blob([script], { type: 'text/plain' });
        this.downloadBlob(blob, `gsea_${plotType}_plot.R`);
    }

    _rEscape(s) { return s.replace(/\\/g, '\\\\').replace(/"/g, '\\"'); }

    /** Wrap a long R vector into ~80-char lines for readability */
    _rWrapVec(items, indent = '        ') {
        let lines = [], line = '';
        for (const item of items) {
            if (line && (line + ', ' + item).length > 80) {
                lines.push(line);
                line = item;
            } else {
                line = line ? line + ', ' + item : item;
            }
        }
        if (line) lines.push(line);
        return lines.join(',\n' + indent);
    }

    _rScriptBubble() {
        const fdrThresh = parseFloat(this.settings.fdrDisplayThreshold);
        const topN = this.settings.topN;
        let filtered = this.results.slice();
        if (fdrThresh < 1) filtered = filtered.filter(r => r.fdr < fdrThresh);
        const top = filtered.sort((a, b) => Math.abs(b.nes) - Math.abs(a.nes)).slice(0, topN);
        top.sort((a, b) => a.nes - b.nes);

        if (top.length === 0) return '# No gene sets pass the current FDR threshold\n';

        const h = Math.max(4, top.length * 0.35 + 1.5).toFixed(1);

        let s = `# Enrichment Overview (Lollipop/Bubble Plot)\n`;
        s += `# Install if needed: install.packages("ggplot2")\n\n`;
        s += `library(ggplot2)\n\n`;

        // Embed data — wrap long vectors for readability
        const nameVec = top.map(r => `"${this._rEscape(this.cleanName(r.name))}"`).join(',\n    ');
        const nesVec  = top.map(r => r.nes.toFixed(4)).join(', ');
        const fdrVec  = top.map(r => r.fdr < 0.001 ? r.fdr.toExponential(4) : r.fdr.toFixed(6)).join(', ');
        const sizeVec = top.map(r => r.size).join(', ');

        s += `df <- data.frame(\n`;
        s += `    name = c(\n    ${nameVec}\n    ),\n`;
        s += `    nes  = c(${nesVec}),\n`;
        s += `    fdr  = c(${fdrVec}),\n`;
        s += `    size = c(${sizeVec}),\n`;
        s += `    stringsAsFactors = FALSE\n)\n\n`;
        s += `# Preserve order from enrichment analysis\n`;
        s += `df$name <- factor(df$name, levels = df$name)\n`;
        s += `df$log_fdr <- -log10(pmax(df$fdr, 1e-10))\n\n`;

        s += `# Plot\n`;
        s += `p <- ggplot(df, aes(x = nes, y = name)) +\n`;
        s += `    geom_segment(aes(x = 0, xend = nes, y = name, yend = name),\n`;
        s += `                 color = "#d1d5db", linewidth = 0.6) +\n`;
        s += `    geom_point(aes(size = size, color = log_fdr)) +\n`;
        s += `    scale_color_gradientn(\n`;
        s += `        colors = c("#ffffb2", "#fd8d3c", "#e31a1c", "#800026"),\n`;
        s += `        name = "FDR\\n(-log10)"\n`;
        s += `    ) +\n`;
        s += `    scale_size_continuous(range = c(3, 10), name = "Gene set\\nsize") +\n`;
        s += `    geom_vline(xintercept = 0, color = "#9ca3af", linewidth = 0.5) +\n`;
        s += `    labs(x = "Normalized Enrichment Score (NES)", y = NULL) +\n`;
        s += `    theme_minimal(base_family = "sans") +\n`;
        s += `    theme(\n`;
        s += `        panel.grid.major.y = element_blank(),\n`;
        s += `        panel.grid.minor   = element_blank(),\n`;
        s += `        axis.text.y        = element_text(size = 9)\n`;
        s += `    )\n\n`;
        s += `print(p)\n`;
        s += `ggsave("gsea_overview.pdf", plot = p, width = 8, height = ${h})\n`;
        s += `cat("Plot saved to:", normalizePath("gsea_overview.pdf"), "\\n")\n`;

        return s;
    }

    _rScriptES() {
        const sel = document.getElementById('geneSetSelector').value;
        const result = this.results.find(r => r.name === sel);
        if (!result) return '# No gene set selected\n';

        const N = this.rankedList.genes.length;
        const metrics = this.rankedList.metrics;
        const cleanName = this.cleanName(result.name);

        // Subsample for R script to avoid huge inline data
        const maxPts = 2000;
        const step = N > maxPts ? Math.ceil(N / maxPts) : 1;
        const ranks = [], esVals = [], metVals = [];
        for (let i = 0; i < N; i += step) {
            ranks.push(i);
            esVals.push(result.runningES[i]);
            metVals.push(metrics[i]);
        }
        if (ranks[ranks.length - 1] !== N - 1) {
            ranks.push(N - 1);
            esVals.push(result.runningES[N - 1]);
            metVals.push(metrics[N - 1]);
        }

        let s = `# Enrichment Score Plot: ${cleanName}\n`;
        s += `# Install if needed: install.packages(c("ggplot2", "patchwork"))\n\n`;
        s += `library(ggplot2)\nlibrary(patchwork)\n\n`;

        // Stats
        s += `# --- Gene set statistics ---\n`;
        s += `gene_set_name <- "${this._rEscape(cleanName)}"\n`;
        s += `nes      <- ${result.nes.toFixed(4)}\n`;
        s += `fdr      <- ${result.fdr < 0.001 ? result.fdr.toExponential(4) : result.fdr.toFixed(6)}\n`;
        s += `pval     <- ${result.pvalue < 0.001 ? result.pvalue.toExponential(4) : result.pvalue.toFixed(6)}\n`;
        s += `set_size <- ${result.size}\n\n`;

        // ES data — wrap long vectors
        s += `# Running enrichment score (subsampled to ${ranks.length} points)\n`;
        s += `es_df <- data.frame(\n`;
        s += `    rank = c(${this._rWrapVec(ranks.map(String))}),\n`;
        s += `    es   = c(${this._rWrapVec(esVals.map(v => v.toFixed(6)))})\n)\n\n`;

        s += `# Hit positions (genes in set)\n`;
        s += `hits <- c(${this._rWrapVec(Array.from(result.hits).map(String))})\n\n`;

        s += `# Ranked list metric (subsampled)\n`;
        s += `met_df <- data.frame(\n`;
        s += `    rank   = c(${this._rWrapVec(ranks.map(String))}),\n`;
        s += `    metric = c(${this._rWrapVec(metVals.map(v => v.toFixed(6)))})\n)\n\n`;

        // ES line color with alpha for fill
        const esCol = this.settings.esLineColor;

        // ES panel
        s += `# --- Panel 1: Enrichment Score curve ---\n`;
        s += `p1 <- ggplot(es_df, aes(x = rank, y = es)) +\n`;
        s += `    geom_area(fill = alpha("${esCol}", 0.12), color = NA) +\n`;
        s += `    geom_line(color = "${esCol}", linewidth = 0.8) +\n`;
        s += `    geom_hline(yintercept = 0, linetype = "dotted", color = "#aaa") +\n`;
        s += `    annotate("label", x = 0, y = Inf, vjust = 1.2, hjust = -0.05,\n`;
        s += `             label = paste0("NES = ", round(nes, 2),\n`;
        s += `                            "\\nFDR = ", formatC(fdr, format = "e", digits = 2),\n`;
        s += `                            "\\np = ", formatC(pval, format = "e", digits = 2),\n`;
        s += `                            "\\nSize = ", set_size),\n`;
        s += `             size = 3, family = "mono", fill = "white", label.size = 0.3) +\n`;
        s += `    labs(y = "Enrichment Score (ES)", title = gene_set_name) +\n`;
        s += `    theme_minimal(base_family = "sans") +\n`;
        s += `    theme(axis.title.x = element_blank(), axis.text.x = element_blank(),\n`;
        s += `          plot.title = element_text(size = 11, face = "bold"),\n`;
        s += `          panel.grid.minor = element_blank())\n\n`;

        // Hit markers panel
        s += `# --- Panel 2: Hit markers ---\n`;
        s += `hit_df <- data.frame(pos = hits)\n`;
        s += `p2 <- ggplot(hit_df, aes(x = pos)) +\n`;
        s += `    geom_segment(aes(x = pos, xend = pos, y = 0, yend = 1), linewidth = 0.3) +\n`;
        s += `    scale_y_continuous(expand = c(0, 0)) +\n`;
        s += `    scale_x_continuous(limits = range(es_df$rank)) +\n`;
        s += `    theme_void() +\n`;
        s += `    theme(plot.margin = margin(0, 5.5, 0, 5.5))\n\n`;

        // Metric panel
        s += `# --- Panel 3: Ranked list metric ---\n`;
        s += `p3 <- ggplot(met_df, aes(x = rank, y = metric, fill = metric > 0)) +\n`;
        s += `    geom_col(width = ${step > 1 ? step * 1.05 : 1}) +\n`;
        s += `    scale_fill_manual(values = c("TRUE" = "${this.settings.positiveColor}",\n`;
        s += `                                  "FALSE" = "${this.settings.negativeColor}"),\n`;
        s += `                      guide = "none") +\n`;
        s += `    geom_hline(yintercept = 0, linewidth = 0.4) +\n`;
        s += `    labs(x = "Rank in Ordered Dataset", y = "Metric") +\n`;
        s += `    theme_minimal(base_family = "sans") +\n`;
        s += `    theme(panel.grid.minor = element_blank())\n\n`;

        // Combine
        const fname = sel.replace(/[^a-zA-Z0-9]/g, '_').substring(0, 40);
        s += `# --- Combine panels ---\n`;
        s += `combined <- p1 / p2 / p3 + plot_layout(heights = c(5, 0.8, 2.5))\n`;
        s += `print(combined)\n`;
        s += `ggsave("gsea_es_${fname}.pdf", plot = combined, width = 7, height = 6)\n`;
        s += `cat("Plot saved to:", normalizePath("gsea_es_${fname}.pdf"), "\\n")\n`;

        return s;
    }

    _rScriptRanked() {
        const N = this.rankedList.genes.length;
        const metrics = this.rankedList.metrics;
        const metricLabel = document.getElementById('metricColumn').value || 'Ranking Metric';

        // Subsample
        const maxBars = 600;
        const step = Math.max(1, Math.ceil(N / maxBars));
        const ranks = [], vals = [];
        for (let i = 0; i < N; i += step) {
            ranks.push(i);
            vals.push(metrics[i]);
        }

        let s = `# Ranked List Metric Plot\n`;
        s += `# Install if needed: install.packages("ggplot2")\n\n`;
        s += `library(ggplot2)\n\n`;

        s += `df <- data.frame(\n`;
        s += `    rank   = c(${this._rWrapVec(ranks.map(String))}),\n`;
        s += `    metric = c(${this._rWrapVec(vals.map(v => v.toFixed(6)))})\n)\n\n`;

        s += `p <- ggplot(df, aes(x = rank, y = metric, fill = metric > 0)) +\n`;
        s += `    geom_col(width = ${step * 1.05}) +\n`;
        s += `    scale_fill_manual(values = c("TRUE" = "${this.settings.positiveColor}",\n`;
        s += `                                  "FALSE" = "${this.settings.negativeColor}"),\n`;
        s += `                      guide = "none") +\n`;
        s += `    geom_hline(yintercept = 0, linewidth = 0.5, color = "#333") +\n`;
        s += `    labs(x = "Rank in Ordered Dataset",\n`;
        s += `         y = "Ranked list metric (${this._rEscape(metricLabel)})") +\n`;
        s += `    theme_minimal(base_family = "sans") +\n`;
        s += `    theme(\n`;
        s += `        panel.grid.minor = element_blank(),\n`;
        s += `        panel.border = element_rect(color = "#333", fill = NA, linewidth = 0.5)\n`;
        s += `    )\n\n`;
        s += `print(p)\n`;
        s += `ggsave("gsea_ranked_list.pdf", plot = p, width = 8, height = 3.5)\n`;
        s += `cat("Plot saved to:", normalizePath("gsea_ranked_list.pdf"), "\\n")\n`;

        return s;
    }

    _rScriptOverlap() {
        // Need to compute the overlap matrix (same as renderOverlapHeatmap)
        const maxSets = parseInt(document.getElementById('overlapMaxSets').value) || 20;
        let sigSets = this.results
            .filter(r => r.fdr < 0.25)
            .sort((a, b) => Math.abs(b.nes) - Math.abs(a.nes))
            .slice(0, maxSets);

        if (sigSets.length < 2) return '# Need at least 2 significant gene sets for overlap\n';

        const activeGeneSets = this.getActiveGeneSets();
        const setGenes = {};
        for (const r of sigSets) {
            if (activeGeneSets[r.name]) {
                const rankedGenesUpper = new Set(this.rankedList.genes.map(g => g.toUpperCase()));
                setGenes[r.name] = activeGeneSets[r.name].map(g => g.toUpperCase()).filter(g => rankedGenesUpper.has(g));
            }
        }
        // Filter to sets we have genes for
        sigSets = sigSets.filter(r => setGenes[r.name] && setGenes[r.name].length > 0);
        if (sigSets.length < 2) return '# Not enough gene sets with gene data for overlap\n';

        const names = sigSets.map(r => this.cleanName(r.name));
        const n = sigSets.length;

        // Compute Jaccard matrix
        const matrix = [];
        for (let i = 0; i < n; i++) {
            const row = [];
            for (let j = 0; j < n; j++) {
                if (i === j) { row.push(1); continue; }
                const setA = new Set(setGenes[sigSets[i].name]);
                const setB = new Set(setGenes[sigSets[j].name]);
                let inter = 0;
                for (const g of setA) if (setB.has(g)) inter++;
                const union = setA.size + setB.size - inter;
                row.push(union > 0 ? +(inter / union).toFixed(4) : 0);
            }
            matrix.push(row);
        }

        const nameVec = names.map(nm => `"${this._rEscape(nm)}"`).join(',\n    ');

        let s = `# Gene Set Overlap Heatmap (Jaccard Similarity)\n`;
        s += `# Install if needed: install.packages("ggplot2")\n\n`;
        s += `library(ggplot2)\n\n`;

        s += `set_names <- c(\n    ${nameVec}\n)\n\n`;

        s += `# Jaccard similarity matrix\n`;
        s += `jaccard <- matrix(c(\n`;
        for (let i = 0; i < n; i++) {
            s += `    ${matrix[i].join(', ')}${i < n - 1 ? ',' : ''}\n`;
        }
        s += `), nrow = ${n}, byrow = TRUE)\n`;
        s += `rownames(jaccard) <- set_names\n`;
        s += `colnames(jaccard) <- set_names\n\n`;

        s += `# Convert to long format (base R — no extra packages needed)\n`;
        s += `df_long <- expand.grid(row = set_names, col = set_names,\n`;
        s += `                       stringsAsFactors = FALSE)\n`;
        s += `df_long$jaccard <- as.vector(t(jaccard))\n`;
        s += `df_long$row <- factor(df_long$row, levels = rev(set_names))\n`;
        s += `df_long$col <- factor(df_long$col, levels = set_names)\n\n`;

        const w = Math.max(6, n * 0.5 + 2).toFixed(1);
        const h = Math.max(5, n * 0.5 + 1.5).toFixed(1);

        s += `p <- ggplot(df_long, aes(x = col, y = row, fill = jaccard)) +\n`;
        s += `    geom_tile(color = "white", linewidth = 0.3) +\n`;
        s += `    scale_fill_gradientn(\n`;
        s += `        colors = c("#ffffff", "#c5e8bc", "#5a9f4a", "#2d5a27"),\n`;
        s += `        name   = "Jaccard\\nIndex",\n`;
        s += `        limits = c(0, 1)\n`;
        s += `    ) +\n`;
        s += `    labs(x = NULL, y = NULL) +\n`;
        s += `    theme_minimal(base_family = "sans") +\n`;
        s += `    theme(\n`;
        s += `        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),\n`;
        s += `        axis.text.y = element_text(size = 8),\n`;
        s += `        panel.grid  = element_blank()\n`;
        s += `    ) +\n`;
        s += `    coord_fixed()\n\n`;
        s += `print(p)\n`;
        s += `ggsave("gsea_overlap.pdf", plot = p, width = ${w}, height = ${h})\n`;
        s += `cat("Plot saved to:", normalizePath("gsea_overlap.pdf"), "\\n")\n`;

        return s;
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
    async loadExampleData(cellLine, dataType) {
        const cellLineInfo = {
            A375: 'Skin cancer (Melanoma)',
            A549: 'Lung cancer (Adenocarcinoma)',
            HT29: 'Colon cancer (Colorectal)',
            Raji: 'Blood cancer (B-cell Lymphoma)',
            U251: 'Brain cancer (Glioblastoma)'
        };
        const metricMap = {
            expression: 'Expression_zscore',
            crispr: 'Chronos_score'
        };
        const labelMap = {
            expression: 'Expression (z-score)',
            crispr: 'CRISPR (Chronos gene effect)'
        };

        this.showStatus('uploadStatus', 'info', `Loading ${cellLine} ${dataType} data from DepMap...`);

        try {
            const resp = await fetch(`web_data/depmap_${dataType}_${cellLine}.json`);
            if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
            const data = await resp.json();

            this.rawData = data;
            const metricCol = metricMap[dataType];
            this.populateColumnDropdowns(['Gene', metricCol]);
            document.getElementById('geneColumn').value = 'Gene';
            document.getElementById('metricColumn').value = metricCol;

            // Auto-set data type (sync both selectors)
            const dtVal = dataType === 'crispr' ? 'crispr' : 'expression';
            document.getElementById('dataType').value = dtVal;
            const dtInline = document.getElementById('dataTypeInline');
            if (dtInline) dtInline.value = dtVal;
            this.settings.dataType = dtVal;

            this.showStatus('uploadStatus', 'success',
                `Loaded ${cellLine} \u2014 ${cellLineInfo[cellLine]} \u2014 ${labelMap[dataType]}: ${data.length.toLocaleString()} genes (DepMap)`);
            this.checkReady();
        } catch (err) {
            this.showStatus('uploadStatus', 'error', `Failed to load example data: ${err.message}`);
        }
    }

    // --------------------------------------------------------
    // Reset Tab Filters
    // --------------------------------------------------------
    resetTabFilters(tab) {
        // Reset hidden sets
        this._hiddenSets.clear();
        this._pinnedBubbleSets.clear();
        this._gsfInitialized = false;
        this._gsfAutoSelectOnOpen = false;

        if (tab === 'overview' || tab === 'all') {
            const el = (id, val) => { const e = document.getElementById(id); if (e) e.value = val; };
            el('overviewFdrFilter', '0.25');
            el('overviewPvalFilter', '1');
            el('overviewDirectionFilter', 'all');
            el('overviewTopN', '20');
            el('overviewOverlapFilter', '0');
            el('fdrDisplayThreshold', '0.25');
            el('pvalDisplayThreshold', '1');
            this._overlapFilterThreshold = 0;
            this.renderBubblePlot();
        }
        if (tab === 'enrichment' || tab === 'all') {
            const el = (id, val) => { const e = document.getElementById(id); if (e) e.value = val; };
            el('esFdrFilter', '0.25');
            el('esPvalFilter', '1');
            el('esDirectionFilter', 'all');
            el('esRedundancyFilter', '0');
            this.populateGeneSetSelector();
            this._initGeneSetSearch();
        }
        if (tab === 'table' || tab === 'all') {
            const el = (id, val) => { const e = document.getElementById(id); if (e) e.value = val; };
            el('fdrFilter', '0.25');
            el('pvalueFilter', 'all');
            el('directionFilter', 'all');
            el('tableOverlapFilter', '0');
            el('tableFilter', '');
            this.filterAndRenderTable();
        }
        if (tab === 'overlap' || tab === 'all') {
            const el = (id, val) => { const e = document.getElementById(id); if (e) e.value = val; };
            el('overlapFdrFilter', '0.25');
            el('overlapPvalFilter', '1');
            el('overlapDirectionFilter', 'all');
            el('overlapRedundancyFilter', '50');
            el('overlapMaxSets', '20');
            this.renderOverlapHeatmap();
        }
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
        this.geneInfoCache = {};
        this._statsBoxCornerCache = {};
        this.selectedGeneSets = new Set();
        this.useCustomSelection = false;
        this._gsbAllItems = [];
        this._gsbFlatList = [];
        this._gsbVisualRows = [];
        this._gsbItemMap = null;

        // Reset overlap filter state
        this._overlapCache = null;
        this._hiddenSets = new Set();
        this._overlapFilterThreshold = 0;

        // Reset UI
        const sidebar = document.querySelector('.sidebar');
        if (sidebar) sidebar.style.paddingTop = '';
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
        document.getElementById('rerunSection').style.display = 'none';
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

        // Reset overlap filter dropdowns
        const bubbleOF = document.getElementById('bubbleOverlapFilter');
        if (bubbleOF) bubbleOF.value = '0';
        const tableOF = document.getElementById('tableOverlapFilter');
        if (tableOF) tableOF.value = '0';

        // Reset checkboxes
        document.getElementById('checkHallmark').checked = true;
        ['checkC2kegg', 'checkC2reactome', 'checkC2wp', 'checkC2biocarta', 'checkC2pid', 'checkC2cgp',
         'checkC5bp', 'checkC5cc', 'checkC5mf', 'checkC5hpo'].forEach(id => {
            document.getElementById(id).checked = false;
        });
        document.getElementById('checkCustomSelection').checked = false;
        document.getElementById('customSelectionCount').textContent = '';

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
    // Gene Set Browser
    // --------------------------------------------------------
    async openGeneSetBrowser() {
        // Auto-load ALL collections (not just checked ones)
        const collections = [
            { id: 'hallmark', file: 'h.all.v2023.2.Hs.json' },
            { id: 'c2kegg', file: 'c2.cp.kegg.v2023.2.Hs.json' },
            { id: 'c2reactome', file: 'c2.cp.reactome.v2023.2.Hs.json' },
            { id: 'c2wp', file: 'c2.cp.wikipathways.v2023.2.Hs.json' },
            { id: 'c2biocarta', file: 'c2.cp.biocarta.v2023.2.Hs.json' },
            { id: 'c2pid', file: 'c2.cp.pid.v2023.2.Hs.json' },
            { id: 'c2cgp', file: 'c2.cgp.v2023.2.Hs.json' },
            { id: 'c3', file: 'c3.all.v2023.2.Hs.json' },
            { id: 'c5bp', file: 'c5.go.bp.v2023.2.Hs.json' },
            { id: 'c5cc', file: 'c5.go.cc.v2023.2.Hs.json' },
            { id: 'c5mf', file: 'c5.go.mf.v2023.2.Hs.json' },
            { id: 'c5hpo', file: 'c5.hpo.v2023.2.Hs.json' },
            { id: 'c6', file: 'c6.all.v2023.2.Hs.json' },
            { id: 'c7', file: 'c7.all.v2023.2.Hs.json' },
            { id: 'c8', file: 'c8.all.v2023.2.Hs.json' }
        ];

        // Show modal immediately with a loading state
        document.getElementById('gsbModal').classList.add('open');
        document.getElementById('gsbBackdrop').classList.add('open');
        document.getElementById('gsbRows').innerHTML = '<div style="padding: 40px 20px; text-align: center; color: #888; font-size: 14px;">Loading all gene set collections...</div>';

        for (const coll of collections) {
            if (!this.geneSets[coll.id]) {
                try {
                    const resp = await fetch('web_data/' + coll.file);
                    if (resp.ok) this.geneSets[coll.id] = await resp.json();
                } catch (e) {
                    console.warn(`Failed to load ${coll.id}:`, e);
                }
            }
        }

        // Build master list
        this._gsbAllItems = [];
        const tierMap = { hallmark: 1, c2: 2, c5: 2, c3: 3, c6: 3, c7: 4, c8: 4, custom: 0 };
        const addCollection = (id, label, data) => {
            if (!data) return;
            for (const [name, genes] of Object.entries(data)) {
                this._gsbAllItems.push({
                    name,
                    collection: id,
                    collectionLabel: label,
                    tier: tierMap[id] || 0,
                    size: genes.length,
                    genes,
                    nameLower: name.toLowerCase(),
                    searchText: name.toLowerCase().replace(/_/g, ' ')
                });
            }
        };
        addCollection('hallmark', 'Hallmark', this.geneSets['hallmark']);
        addCollection('c2kegg', 'C2: KEGG', this.geneSets['c2kegg']);
        addCollection('c2reactome', 'C2: Reactome', this.geneSets['c2reactome']);
        addCollection('c2wp', 'C2: WikiPathways', this.geneSets['c2wp']);
        addCollection('c2biocarta', 'C2: BioCarta', this.geneSets['c2biocarta']);
        addCollection('c2pid', 'C2: PID', this.geneSets['c2pid']);
        addCollection('c2cgp', 'C2: CGP', this.geneSets['c2cgp']);
        addCollection('c3', 'C3: Regulatory', this.geneSets['c3']);
        addCollection('c5bp', 'C5: GO:BP', this.geneSets['c5bp']);
        addCollection('c5cc', 'C5: GO:CC', this.geneSets['c5cc']);
        addCollection('c5mf', 'C5: GO:MF', this.geneSets['c5mf']);
        addCollection('c5hpo', 'C5: HPO', this.geneSets['c5hpo']);
        addCollection('c6', 'C6: Oncogenic', this.geneSets['c6']);
        addCollection('c7', 'C7: Immunologic', this.geneSets['c7']);
        addCollection('c8', 'C8: Cell Type', this.geneSets['c8']);
        if (this.customGeneSets && Object.keys(this.customGeneSets).length > 0) {
            addCollection('custom', 'Custom', this.customGeneSets);
        }

        // Sort: by collection order, then alphabetically within each
        const collOrder = { hallmark: 0, c2kegg: 1, c2reactome: 2, c2wp: 3, c2biocarta: 4, c2pid: 5, c2cgp: 6, c3: 7, c5bp: 8, c5cc: 9, c5mf: 10, c5hpo: 11, c6: 12, c7: 13, c8: 14, custom: 15 };
        this._gsbAllItems.sort((a, b) => {
            if (collOrder[a.collection] !== collOrder[b.collection]) {
                return collOrder[a.collection] - collOrder[b.collection];
            }
            return a.nameLower.localeCompare(b.nameLower);
        });

        // Build lookup map for hover detail (O(1) instead of linear search)
        this._gsbItemMap = new Map();
        for (const item of this._gsbAllItems) {
            this._gsbItemMap.set(item.name, item);
        }

        // Build gene → set names index for gene-based search
        this._gsbGeneIndex = new Map();
        for (const item of this._gsbAllItems) {
            for (const gene of item.genes) {
                const g = gene.toUpperCase();
                if (!this._gsbGeneIndex.has(g)) this._gsbGeneIndex.set(g, []);
                this._gsbGeneIndex.get(g).push(item.name);
            }
        }

        // Pre-compute gene Sets and diversity caches for instant Jaccard filtering
        // (runs synchronously — typically <500ms even for 30k+ sets)
        this._gsbPrecomputeDiversity();

        // If first time opening (no prior custom selection), start with nothing selected
        if (!this.useCustomSelection) {
            this.selectedGeneSets = new Set();
        }

        // Reset search and filter
        document.getElementById('gsbSearch').value = '';
        document.getElementById('gsbTierFilter').value = 'all';
        this._updateCollectionFilterForTier();
        document.getElementById('gsbCollectionFilter').value = 'all';
        document.getElementById('gsbJaccardFilter').value = 'off';

        // Reset detail panel
        document.getElementById('gsbDetailTitle').textContent = 'Gene Set Details';
        document.getElementById('gsbDetailBody').innerHTML = '<p class="gsb-detail-hint">Hover over any gene set to see its details here.</p>';

        // Apply filter and render
        this.filterGeneSetBrowser();
        this.updateGsbSelectionCount();
        document.getElementById('gsbSearch').focus();
    }

    closeGeneSetBrowser() {
        document.getElementById('gsbModal').classList.remove('open');
        document.getElementById('gsbBackdrop').classList.remove('open');
    }

    _gsbPrecomputeDiversity() {
        const thresholds = [0.1, 0.2, 0.3, 0.5];
        // Build uppercase gene Sets for all items
        for (const item of this._gsbAllItems) {
            if (!item._geneSet) item._geneSet = new Set(item.genes.map(g => g.toUpperCase()));
        }
        // Group items by collection
        const byCollection = {};
        for (const item of this._gsbAllItems) {
            if (!byCollection[item.collection]) byCollection[item.collection] = [];
            byCollection[item.collection].push(item);
        }
        // For each collection, run greedy diversity selection at each threshold
        // Cache: _gsbDiversityCache[threshold] = Set of item names that pass
        this._gsbDiversityCache = {};
        for (const t of thresholds) {
            const passNames = new Set();
            for (const items of Object.values(byCollection)) {
                const keptGeneSets = [];
                for (const item of items) {
                    const itemGenes = item._geneSet;
                    let tooSimilar = false;
                    for (const keptGenes of keptGeneSets) {
                        const sA = itemGenes.size, sB = keptGenes.size;
                        const minSize = Math.min(sA, sB), maxSize = Math.max(sA, sB);
                        if (minSize / maxSize < t) continue;
                        let intersection = 0;
                        const smaller = sA <= sB ? itemGenes : keptGenes;
                        const larger = sA <= sB ? keptGenes : itemGenes;
                        for (const g of smaller) {
                            if (larger.has(g)) intersection++;
                        }
                        const union = sA + sB - intersection;
                        if (union > 0 && intersection / union >= t) {
                            tooSimilar = true;
                            break;
                        }
                    }
                    if (!tooSimilar) {
                        passNames.add(item.name);
                        keptGeneSets.push(itemGenes);
                    }
                }
            }
            this._gsbDiversityCache[t] = passNames;
        }
    }

    _updateCollectionFilterForTier() {
        const tierFilter = document.getElementById('gsbTierFilter').value;
        const collEl = document.getElementById('gsbCollectionFilter');
        const tierCollections = {
            'all': null,
            '1': ['hallmark'],
            '2': ['c2kegg', 'c2reactome', 'c2wp', 'c2biocarta', 'c2pid', 'c2cgp', 'c5bp', 'c5cc', 'c5mf', 'c5hpo'],
            '3': ['c3', 'c6'],
            '4': ['c7', 'c8']
        };
        const allOptions = [
            { value: 'all', label: 'All collections' },
            { value: 'hallmark', label: 'H: Hallmark' },
            { value: 'c2kegg', label: 'C2: KEGG' },
            { value: 'c2reactome', label: 'C2: Reactome' },
            { value: 'c2wp', label: 'C2: WikiPathways' },
            { value: 'c2biocarta', label: 'C2: BioCarta' },
            { value: 'c2pid', label: 'C2: PID' },
            { value: 'c2cgp', label: 'C2: CGP' },
            { value: 'c3', label: 'C3: Regulatory' },
            { value: 'c5bp', label: 'C5: GO:BP' },
            { value: 'c5cc', label: 'C5: GO:CC' },
            { value: 'c5mf', label: 'C5: GO:MF' },
            { value: 'c5hpo', label: 'C5: HPO' },
            { value: 'c6', label: 'C6: Oncogenic' },
            { value: 'c7', label: 'C7: Immunologic' },
            { value: 'c8', label: 'C8: Cell Type' },
            { value: 'custom', label: 'Custom' }
        ];
        const allowed = tierCollections[tierFilter];
        collEl.innerHTML = '';
        for (const opt of allOptions) {
            if (opt.value === 'all' || !allowed || allowed.includes(opt.value)) {
                const o = document.createElement('option');
                o.value = opt.value;
                o.textContent = opt.label;
                collEl.appendChild(o);
            }
        }
        collEl.value = 'all';
    }

    filterGeneSetBrowser() {
        const query = document.getElementById('gsbSearch').value.toLowerCase().trim();
        const collFilter = document.getElementById('gsbCollectionFilter').value;
        const tierFilter = document.getElementById('gsbTierFilter').value;

        // Filter items
        let filtered = this._gsbAllItems;
        if (tierFilter !== 'all') {
            const tierCollections = { '1': ['hallmark'], '2': ['c2kegg', 'c2reactome', 'c2wp', 'c2biocarta', 'c2pid', 'c2cgp', 'c5bp', 'c5cc', 'c5mf', 'c5hpo'], '3': ['c3', 'c6'], '4': ['c7', 'c8'] };
            const allowed = new Set(tierCollections[tierFilter] || []);
            filtered = filtered.filter(item => allowed.has(item.collection));
        }
        if (collFilter !== 'all') {
            filtered = filtered.filter(item => item.collection === collFilter);
        }
        if (query) {
            // 1) Name search: match full query against cleaned gene set names
            const nameMatchSet = new Set();
            for (const item of filtered) {
                if (item.searchText.includes(query)) nameMatchSet.add(item.name);
            }

            // 2) Gene search: split query by commas/spaces into tokens,
            //    find all gene sets containing any of those genes
            const geneMatchSet = new Set();
            const tokens = query.split(/[\s,]+/).map(t => t.trim().toUpperCase()).filter(t => t.length >= 2);
            for (const token of tokens) {
                const sets = this._gsbGeneIndex.get(token);
                if (sets) for (const name of sets) geneMatchSet.add(name);
            }

            // Union both result sets
            filtered = filtered.filter(item =>
                nameMatchSet.has(item.name) || geneMatchSet.has(item.name)
            );

            // Store for display hint
            this._gsbLastGeneHits = geneMatchSet.size > 0 ? tokens.filter(t => this._gsbGeneIndex.has(t)) : [];
        } else {
            this._gsbLastGeneHits = [];
        }

        // Jaccard diversity filter: use pre-computed cache for instant filtering
        const jaccardThreshold = document.getElementById('gsbJaccardFilter').value;
        let jaccardFiltered = 0;
        if (jaccardThreshold !== 'off' && this._gsbDiversityCache) {
            const maxJ = parseFloat(jaccardThreshold);
            const passNames = this._gsbDiversityCache[maxJ];
            if (passNames) {
                const before = filtered.length;
                filtered = filtered.filter(item => passNames.has(item.name));
                jaccardFiltered = before - filtered.length;
            }
        }

        // Build flat list with section headers
        this._gsbFlatList = [];

        // If gene matches found, show info banner
        if (this._gsbLastGeneHits.length > 0) {
            this._gsbFlatList.push({
                type: 'info',
                text: `Gene match: ${this._gsbLastGeneHits.join(', ')} found in ${filtered.length} gene sets`
            });
        }

        // If Jaccard filter active, show info banner
        if (jaccardFiltered > 0) {
            this._gsbFlatList.push({
                type: 'info',
                text: `Diversity filter: showing ${filtered.length} diverse sets (${jaccardFiltered} similar sets hidden, J ≥ ${jaccardThreshold})`
            });
        }

        let lastCollection = null;
        const collectionCounts = {};

        // Pre-count per collection
        for (const item of filtered) {
            collectionCounts[item.collection] = (collectionCounts[item.collection] || 0) + 1;
        }

        for (const item of filtered) {
            if (item.collection !== lastCollection) {
                this._gsbFlatList.push({
                    type: 'header',
                    collection: item.collection,
                    label: item.collectionLabel,
                    tier: item.tier,
                    count: collectionCounts[item.collection]
                });
                lastCollection = item.collection;
            }
            this._gsbFlatList.push({ type: 'item', ...item });
        }

        // Build visual rows: headers span full width, items grouped into N columns
        const cols = this._gsbColumns;
        this._gsbVisualRows = [];
        let pendingItems = [];

        const flushItems = () => {
            while (pendingItems.length > 0) {
                this._gsbVisualRows.push({
                    type: 'items',
                    items: pendingItems.splice(0, cols)
                });
            }
        };

        for (const entry of this._gsbFlatList) {
            if (entry.type === 'header' || entry.type === 'info') {
                flushItems();
                this._gsbVisualRows.push(entry);
            } else {
                pendingItems.push(entry);
            }
        }
        flushItems();

        // Update results count
        const totalItems = this._gsbFlatList.filter(e => e.type === 'item').length;
        const countEl = document.getElementById('gsbResultsCount');
        if (query || collFilter !== 'all') {
            countEl.textContent = `${totalItems.toLocaleString()} gene sets found`;
        } else {
            countEl.textContent = `${totalItems.toLocaleString()} gene sets`;
        }

        // Update scroll spacer height
        const container = document.getElementById('gsbListContainer');
        const spacer = document.getElementById('gsbScrollSpacer');
        spacer.style.height = (this._gsbVisualRows.length * this._gsbRowHeight) + 'px';

        // Reset scroll
        container.scrollTop = 0;
        this.gsbRenderVisibleRows();
    }

    gsbRenderVisibleRows() {
        const container = document.getElementById('gsbListContainer');
        const rowsEl = document.getElementById('gsbRows');
        const scrollTop = container.scrollTop;
        const viewportHeight = container.clientHeight;
        const rowH = this._gsbRowHeight;
        const buffer = this._gsbBufferRows;
        const cols = this._gsbColumns;
        const colPct = (100 / cols).toFixed(4);

        const startIndex = Math.max(0, Math.floor(scrollTop / rowH) - buffer);
        const endIndex = Math.min(
            this._gsbVisualRows.length,
            Math.ceil((scrollTop + viewportHeight) / rowH) + buffer
        );

        let html = '';
        for (let i = startIndex; i < endIndex; i++) {
            const vrow = this._gsbVisualRows[i];
            const top = i * rowH;

            if (vrow.type === 'info') {
                html += `<div class="gsb-info-banner" style="position:absolute; top:${top}px; left:0; right:0; height:${rowH}px;">
                    ${vrow.text}
                </div>`;
            } else if (vrow.type === 'header') {
                const tierNames = { 1: 'Core', 2: 'Standard', 3: 'Specialized', 4: 'Extended' };
                const tierBadge = vrow.tier ? `<span class="gsb-tier-badge tier-${vrow.tier}">${tierNames[vrow.tier]}</span>` : '';
                html += `<div class="gsb-section-header" style="position:absolute; top:${top}px; left:0; right:0; height:${rowH}px;">
                    ${vrow.label} (${vrow.count.toLocaleString()} gene sets)${tierBadge}
                </div>`;
            } else {
                // Render N items side by side
                for (let c = 0; c < vrow.items.length; c++) {
                    const item = vrow.items[c];
                    const checked = this.selectedGeneSets.has(item.name) ? 'checked' : '';
                    const selectedClass = this.selectedGeneSets.has(item.name) ? ' selected' : '';
                    const escapedName = this._escapeAttr(item.name);
                    const leftPct = (c * 100 / cols).toFixed(4);
                    html += `<div class="gsb-row${selectedClass}" data-name="${escapedName}" style="position:absolute; top:${top}px; left:${leftPct}%; width:${colPct}%; height:${rowH}px;">
                        <input type="checkbox" ${checked} data-name="${escapedName}">
                        <span class="gsb-row-name" title="${escapedName}">${this.cleanName(item.name)}</span>
                        <span class="gsb-row-size">${item.size}</span>
                    </div>`;
                }
            }
        }

        rowsEl.innerHTML = html;
    }

    gsbShowDetail(item) {
        document.getElementById('gsbDetailTitle').textContent = this.cleanName(item.name);

        const genes = item.genes;
        const msigdbUrl = `https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/${encodeURIComponent(item.name)}.html`;

        let html = `
            <div class="gsb-detail-stat">
                <span class="gsb-detail-stat-label">Full name</span>
            </div>
            <div style="font-family: 'Roboto Mono', monospace; font-size: 11px; color: #374151; margin-bottom: 10px; word-break: break-all;">
                ${this._escapeAttr(item.name)}
            </div>
            <div class="gsb-detail-stat">
                <span class="gsb-detail-stat-label">Collection</span>
                <span class="gsb-detail-stat-value">${item.collectionLabel}</span>
            </div>
            <div class="gsb-detail-stat">
                <span class="gsb-detail-stat-label">Size</span>
                <span class="gsb-detail-stat-value">${item.size} genes</span>
            </div>
            <div style="margin-top: 8px;">
                <a href="${msigdbUrl}" target="_blank" rel="noopener" style="display: inline-flex; align-items: center; gap: 4px; font-size: 12px; color: var(--green-600); text-decoration: none; font-weight: 600; padding: 4px 8px; background: var(--green-50); border: 1px solid var(--green-200); border-radius: 4px;">🔗 View on MSigDB</a>
            </div>
            <div style="margin-top: 12px;">
                <strong style="font-size: 12px; color: #374151;">Genes (${genes.length}):</strong>
            </div>
            <div class="gsb-detail-genes">
                ${genes.join(', ')}
            </div>
        `;

        document.getElementById('gsbDetailBody').innerHTML = html;
    }

    gsbSelectAllVisible() {
        for (const item of this._gsbFlatList) {
            if (item.type === 'item') {
                this.selectedGeneSets.add(item.name);
            }
        }
        this.updateGsbSelectionCount();
        this.gsbRenderVisibleRows();
    }

    gsbDeselectAll() {
        this.selectedGeneSets.clear();
        this.updateGsbSelectionCount();
        this.gsbRenderVisibleRows();
    }

    updateGsbSelectionCount() {
        const count = this.selectedGeneSets.size;
        document.getElementById('gsbSelectionCount').textContent =
            `${count.toLocaleString()} selected`;

        // Show/hide warning banner for large selections
        const warningEl = document.getElementById('gsbWarningBanner');
        if (warningEl) {
            if (count > 1000) {
                warningEl.style.display = '';
                warningEl.textContent = `⚠ You have selected ${count.toLocaleString()} gene sets. ` +
                    `Analysis with many sets takes longer and produces stricter FDR corrections. ` +
                    `Consider starting with Hallmark (50 sets) for a quick overview.`;
            } else {
                warningEl.style.display = 'none';
            }
        }
    }

    applyGeneSetSelection() {
        if (this.selectedGeneSets.size === 0) {
            alert('No gene sets selected. Please select at least one gene set.');
            return;
        }
        this.useCustomSelection = true;

        // Uncheck collection checkboxes — the browser selection takes over
        ['checkHallmark', 'checkC2kegg', 'checkC2reactome', 'checkC2wp', 'checkC2biocarta', 'checkC2pid', 'checkC2cgp', 'checkC3', 'checkC5bp', 'checkC5cc', 'checkC5mf', 'checkC5hpo', 'checkC6', 'checkC7', 'checkC8'].forEach(id => {
            document.getElementById(id).checked = false;
        });

        // Check the Custom selection row and show count
        document.getElementById('checkCustomSelection').checked = true;
        document.getElementById('customSelectionCount').textContent = this.selectedGeneSets.size.toLocaleString() + ' sets';

        this.closeGeneSetBrowser();
        this.updateGeneSetStatus();
        this.checkReady();
    }

    gsbDownloadSelectedCSV() {
        if (this.selectedGeneSets.size === 0) {
            alert('No gene sets selected.');
            return;
        }

        const rows = [['Gene_Set_Name', 'Collection', 'Size', 'Genes']];

        for (const item of this._gsbAllItems) {
            if (!this.selectedGeneSets.has(item.name)) continue;
            rows.push([
                item.name,
                item.collectionLabel,
                item.size,
                item.genes.join('; ')
            ]);
        }

        const csv = rows.map(row =>
            row.map(val => {
                const s = String(val);
                return (s.includes(',') || s.includes(';') || s.includes('"'))
                    ? '"' + s.replace(/"/g, '""') + '"'
                    : s;
            }).join(',')
        ).join('\n');

        this.downloadBlob(new Blob([csv], { type: 'text/csv' }), 'selected_gene_sets.csv');
    }

    _escapeAttr(str) {
        return str.replace(/&/g, '&amp;').replace(/"/g, '&quot;').replace(/'/g, '&#39;').replace(/</g, '&lt;');
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
        this.settings.showESIndicator = cb('showESIndicator');
        this.settings.showCorrelationLabels = cb('showCorrelationLabels');
        this.settings.showPanelBorders = cb('showPanelBorders');
        this.settings.showHitMarkers = cb('showHitMarkers');
        this.settings.showMetricFill = cb('showMetricFill');
        // Per-graph color settings
        const posColorEl = document.getElementById('positiveColor');
        if (posColorEl) this.settings.positiveColor = posColorEl.value;
        const negColorEl = document.getElementById('negativeColor');
        if (negColorEl) this.settings.negativeColor = negColorEl.value;
        const overlapSchemeEl = document.getElementById('overlapColorScheme');
        if (overlapSchemeEl) this.settings.overlapColorScheme = overlapSchemeEl.value;
        const dataTypeEl = document.getElementById('dataType');
        if (dataTypeEl) this.settings.dataType = dataTypeEl.value;
        // ES panel proportions
        const esPanelESEl = document.getElementById('esPanelES');
        if (esPanelESEl) {
            this.settings.esPanelES = parseInt(esPanelESEl.value) || 55;
            document.getElementById('esPanelESLabel').textContent = this.settings.esPanelES + '%';
        }
        const esPanelHitsEl = document.getElementById('esPanelHits');
        if (esPanelHitsEl) {
            this.settings.esPanelHits = parseInt(esPanelHitsEl.value) || 8;
            document.getElementById('esPanelHitsLabel').textContent = this.settings.esPanelHits + '%';
        }
        const esPanelMetricEl = document.getElementById('esPanelMetric');
        if (esPanelMetricEl) {
            this.settings.esPanelMetric = parseInt(esPanelMetricEl.value) || 28;
            document.getElementById('esPanelMetricLabel').textContent = this.settings.esPanelMetric + '%';
        }
        // Plot heights
        const esHeightEl = document.getElementById('esPlotHeight');
        if (esHeightEl) this.settings.esPlotHeight = parseInt(esHeightEl.value) || 580;
    }

    updateSettings() {
        this.readSettings();
        // Sync gear popup values to overview toolbar
        const ovFdr = document.getElementById('overviewFdrFilter');
        const ovPval = document.getElementById('overviewPvalFilter');
        const ovTopN = document.getElementById('overviewTopN');
        if (ovFdr) ovFdr.value = this.settings.fdrDisplayThreshold;
        if (ovPval) ovPval.value = document.getElementById('pvalDisplayThreshold')?.value || '1';
        if (ovTopN) ovTopN.value = this.settings.topN;
        if (this.results) {
            // Re-render active tab's plots
            if (this.activeTab === 'overview') {
                this.renderBubblePlot();
                this.renderRankedPlot();
            } else if (this.activeTab === 'enrichment') {
                const selected = document.getElementById('geneSetSelector').value;
                if (selected) this.renderESPlot(selected);
            } else if (this.activeTab === 'overlap') {
                this.renderOverlapHeatmap();
            }
            // Also re-render any open tab's plots if global settings changed (font, etc.)
        }
    }

    updateSettingsTabVisibility() {
        // No-op — settings are now per-graph via inline popups
    }

    // --------------------------------------------------------
    // Per-Text-Element Settings
    // --------------------------------------------------------
    _initTextFonts() {
        this.textFonts = {};
        for (const [plotType, elements] of Object.entries(this.textElements)) {
            for (const el of elements) {
                this.textFonts[`${plotType}_${el.key}`] = {
                    family: el.defaultFamily || 'Open Sans',
                    size: el.defaultSize || 12,
                    bold: false,
                    italic: false,
                    visible: true,
                    text: el.defaultText || ''
                };
            }
        }
    }

    _getTextFont(plotType, key) {
        const id = `${plotType}_${key}`;
        const f = this.textFonts[id];
        if (!f) return { family: this.settings.fontFamily + ', sans-serif', size: this.settings.fontSize, wrap: t => t, visible: true, text: '' };
        return {
            family: f.family + ', sans-serif',
            size: f.size,
            visible: f.visible,
            text: f.text,
            wrap: (t) => {
                let r = t;
                if (f.bold) r = `<b>${r}</b>`;
                if (f.italic) r = `<i>${r}</i>`;
                return r;
            }
        };
    }

    openTextSettings(plotType) {
        const panel = document.getElementById('textSettingsPanel');
        const body = document.getElementById('textSettingsPanelBody');
        const elements = this.textElements[plotType];
        if (!elements) return;

        this._activeTextPlotType = plotType;
        const plotLabels = { bubble: 'Overview Plot', ranked: 'Ranked Plot', es: 'Enrichment Plot', overlap: 'Overlap Heatmap' };

        const fontOptions = ['Open Sans', 'Arial', 'Helvetica', 'Times New Roman', 'Georgia', 'Roboto Mono']
            .map(f => `<option value="${f}">${f}</option>`).join('');

        // Check if all elements are bold for the toggle state
        const allBold = elements.every(el => {
            const f = this.textFonts[`${plotType}_${el.key}`];
            return f && f.bold;
        });
        const scaleOffset = this._textScaleOffset[plotType] || 0;
        const scaleDisplay = scaleOffset > 0 ? `+${scaleOffset}` : scaleOffset === 0 ? '0' : `${scaleOffset}`;

        let html = `<div style="font-size: 11px; color: var(--gray-500); margin-bottom: 8px; font-weight: 600;">${plotLabels[plotType] || plotType}</div>`;

        // Scale all controls
        html += `<div class="ts-row" style="background: var(--green-50); border-radius: 4px; padding: 4px 6px; margin-bottom: 8px; border: 1px solid var(--green-200);">`;
        html += `<div class="ts-row-label"><span class="ts-label" style="font-weight: 600;">Scale all</span></div>`;
        html += `<div class="ts-controls">`;
        html += `<button class="ts-btn" data-action="scale-down" data-id="_all" title="Decrease all sizes">&minus;</button>`;
        html += `<span style="display: inline-block; width: 28px; text-align: center; font-size: 11px; font-family: Roboto Mono, monospace;" id="tsScaleValue">${scaleDisplay}</span>`;
        html += `<button class="ts-btn" data-action="scale-up" data-id="_all" title="Increase all sizes">+</button>`;
        html += `<button class="ts-btn ${allBold ? 'ts-btn-active' : ''}" data-action="bold-all" data-id="_all" style="font-weight:bold;" title="Toggle bold on all elements">B</button>`;
        html += `</div></div>`;

        for (const el of elements) {
            const id = `${plotType}_${el.key}`;
            const f = this.textFonts[id];
            if (!f) continue;

            html += `<div class="ts-row" data-ts-id="${id}">`;
            html += `<div class="ts-row-label">`;
            html += `<span class="ts-vis ${f.visible ? '' : 'ts-vis-off'}" data-action="visibility" data-id="${id}" title="Show/Hide">&#128065;</span>`;
            html += `<span class="ts-label">${el.label}</span>`;
            html += `</div>`;

            // Text input for editable elements
            if (el.editable) {
                html += `<input type="text" class="ts-text" data-action="text" data-id="${id}" value="${this._escapeAttr(f.text)}">`;
            }

            // Font controls row
            html += `<div class="ts-controls">`;
            html += `<select class="ts-font-select" data-action="family" data-id="${id}">${fontOptions.replace(`value="${f.family}"`, `value="${f.family}" selected`)}</select>`;
            html += `<input type="number" class="ts-size" data-action="size" data-id="${id}" value="${f.size}" min="6" max="30">`;
            html += `<button class="ts-btn ${f.bold ? 'ts-btn-active' : ''}" data-action="bold" data-id="${id}" style="font-weight:bold;">B</button>`;
            html += `<button class="ts-btn ${f.italic ? 'ts-btn-active' : ''}" data-action="italic" data-id="${id}" style="font-style:italic;">I</button>`;
            html += `</div></div>`;
        }

        body.innerHTML = html;
        panel.style.display = 'block';

        // Event delegation
        if (!this._tsEventsBound) {
            this._tsEventsBound = true;
            body.addEventListener('input', (e) => this._handleTextSettingsInput(e));
            body.addEventListener('change', (e) => this._handleTextSettingsInput(e));
            body.addEventListener('click', (e) => this._handleTextSettingsClick(e));
            // Drag
            this._initPanelDrag('textSettingsPanel', 'textSettingsDragHandle');
        }
    }

    _handleTextSettingsInput(e) {
        const action = e.target.dataset.action;
        const id = e.target.dataset.id;
        if (!action || !id || !this.textFonts[id]) return;
        if (action === 'text') this.textFonts[id].text = e.target.value;
        else if (action === 'family') this.textFonts[id].family = e.target.value;
        else if (action === 'size') this.textFonts[id].size = parseInt(e.target.value) || 12;
        this.updateSettings();
    }

    _handleTextSettingsClick(e) {
        const btn = e.target.closest('[data-action]');
        if (!btn) return;
        const action = btn.dataset.action;
        const id = btn.dataset.id;

        // Prevent document click handler from closing the panel
        // (openTextSettings rebuilds innerHTML, detaching e.target from DOM,
        //  so .closest('.text-settings-panel') returns null on the bubble phase)
        e.stopPropagation();

        // Scale all / bold all actions
        if (action === 'scale-up' || action === 'scale-down' || action === 'bold-all') {
            const pt = this._activeTextPlotType;
            if (!pt || !this.textElements[pt]) return;

            if (action === 'scale-up' || action === 'scale-down') {
                const delta = action === 'scale-up' ? 1 : -1;
                this._textScaleOffset[pt] = (this._textScaleOffset[pt] || 0) + delta;
                for (const el of this.textElements[pt]) {
                    const fid = `${pt}_${el.key}`;
                    if (this.textFonts[fid]) {
                        this.textFonts[fid].size = Math.max(6, Math.min(30, this.textFonts[fid].size + delta));
                    }
                }
            } else if (action === 'bold-all') {
                const allBold = this.textElements[pt].every(el => {
                    const f = this.textFonts[`${pt}_${el.key}`];
                    return f && f.bold;
                });
                const newBold = !allBold;
                for (const el of this.textElements[pt]) {
                    const fid = `${pt}_${el.key}`;
                    if (this.textFonts[fid]) this.textFonts[fid].bold = newBold;
                }
            }

            // Refresh the panel to reflect changes
            this.openTextSettings(pt);
            this.updateSettings();
            return;
        }

        if (!id || !this.textFonts[id]) return;
        if (action === 'bold') {
            this.textFonts[id].bold = !this.textFonts[id].bold;
            btn.classList.toggle('ts-btn-active');
        } else if (action === 'italic') {
            this.textFonts[id].italic = !this.textFonts[id].italic;
            btn.classList.toggle('ts-btn-active');
        } else if (action === 'visibility') {
            this.textFonts[id].visible = !this.textFonts[id].visible;
            btn.classList.toggle('ts-vis-off');
        }
        this.updateSettings();
    }

    closeTextSettings() {
        document.getElementById('textSettingsPanel').style.display = 'none';
    }

    _initPanelDrag(panelId, handleId) {
        const panel = document.getElementById(panelId);
        const handle = document.getElementById(handleId);
        if (!panel || !handle) return;
        let isDragging = false, startX, startY, startLeft, startTop;
        handle.addEventListener('mousedown', (e) => {
            if (e.target.tagName === 'BUTTON') return;
            isDragging = true;
            const rect = panel.getBoundingClientRect();
            startX = e.clientX; startY = e.clientY;
            startLeft = rect.left; startTop = rect.top;
            e.preventDefault();
        });
        document.addEventListener('mousemove', (e) => {
            if (!isDragging) return;
            panel.style.left = (startLeft + (e.clientX - startX)) + 'px';
            panel.style.top = (startTop + (e.clientY - startY)) + 'px';
            panel.style.right = 'auto';
        });
        document.addEventListener('mouseup', () => { isDragging = false; });
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

    toggleGraphSettings(panelId) {
        const panel = document.getElementById(panelId);
        if (!panel) return;
        // Close other open popups first
        document.querySelectorAll('.graph-settings-popup.open').forEach(p => {
            if (p.id !== panelId) p.classList.remove('open');
        });
        panel.classList.toggle('open');
    }

    getPlotDimensions(plotId, defaultW, defaultH) {
        // Per-graph dimension inputs — try graph-specific first, then fallback
        let wId, hId;
        if (plotId === 'bubblePlot') { wId = 'plotWidth'; hId = 'plotHeight'; }
        else if (plotId === 'rankedPlot') { wId = 'rankedPlotWidth'; hId = 'rankedPlotHeight'; }
        else if (plotId === 'esPlot') { wId = 'esPlotWidth'; hId = 'esPlotHeight'; }
        else if (plotId === 'overlapHeatmap') { wId = 'overlapPlotWidth'; hId = 'overlapPlotHeight'; }
        const wEl = wId ? document.getElementById(wId) : null;
        const hEl = hId ? document.getElementById(hId) : null;
        const w = wEl ? (parseInt(wEl.value) || 0) : 0;
        const h = hEl ? (parseInt(hEl.value) || 0) : 0;
        return {
            width: w > 0 ? w : undefined,
            height: h > 0 ? h : defaultH
        };
    }

    /** Sync sidebar top padding to align with main panel card headers (below tab bar) */
    _syncSidebarOffset() {
        requestAnimationFrame(() => {
            const tabBar = document.getElementById('tabBar');
            const sidebar = document.querySelector('.sidebar');
            if (tabBar && sidebar) {
                const h = tabBar.offsetHeight + 8; // 8 = tab-bar margin-bottom
                sidebar.style.paddingTop = h + 'px';
            }
        });
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

    _hexToRgba(hex, alpha) {
        const r = parseInt(hex.slice(1, 3), 16);
        const g = parseInt(hex.slice(3, 5), 16);
        const b = parseInt(hex.slice(5, 7), 16);
        return `rgba(${r},${g},${b},${alpha})`;
    }

    _parseHex(hex) {
        return [
            parseInt(hex.slice(1, 3), 16),
            parseInt(hex.slice(3, 5), 16),
            parseInt(hex.slice(5, 7), 16)
        ];
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
