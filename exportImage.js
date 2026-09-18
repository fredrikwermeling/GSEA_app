//
// Enrich, image export. Ported from Correlate so the two apps export figures
// the same way: one "Image..." button per figure opens a dialog that asks for
// format (PNG, SVG, PDF, TIFF, PowerPoint), print size in cm, DPI and
// background. The figure is taken from Plotly as SVG at its on-screen size,
// so it looks exactly like the screen, then written at the requested physical
// size: PNG/TIFF carry the DPI, PDF and PowerPoint are vector.
// MIT Open source
//

function exportStamp() {
    const d = new Date();
    const p2 = (n) => String(n).padStart(2, '0');
    return `${d.getFullYear()}-${p2(d.getMonth() + 1)}-${p2(d.getDate())}_${p2(d.getHours())}${p2(d.getMinutes())}`;
}

Object.assign(GSEAApp.prototype, {

    // Entry point from the card headers. plotId is the Plotly div.
    async exportImage(plotId) {
        const plotEl = document.getElementById(plotId);
        if (!plotEl || !plotEl.data) { alert('Nothing to export yet. Run an analysis first.'); return; }
        const fl = plotEl._fullLayout || {};
        const w = fl.width || plotEl.offsetWidth || 800;
        const h = fl.height || plotEl.offsetHeight || 500;
        const stems = { bubblePlot: 'enrich_bubble_plot', rankedPlot: 'enrich_ranked_list', esPlot: 'enrich_enrichment_plot', overlapHeatmap: 'enrich_overlap', gsHeatmap: 'enrich_geneset_heatmap' };
        let filename = stems[plotId] || `enrich_${plotId}`;
        if (plotId === 'esPlot') {
            const sel = document.getElementById('geneSetSelector');
            if (sel && sel.value) filename += '_' + this.cleanName(sel.value).replace(/[^A-Za-z0-9]+/g, '_').replace(/^_|_$/g, '').slice(0, 60);
        }
        const dlg = await this._showExportDialog({ plotW: w, plotH: h });
        if (!dlg) return;
        const { background } = dlg;

        const svgDataUrl = await Plotly.toImage(plotEl, { format: 'svg', width: w, height: h });
        let svgStr = svgDataUrl.indexOf('base64,') > -1
            ? atob(svgDataUrl.split('base64,')[1])
            : decodeURIComponent(svgDataUrl.split(',').slice(1).join(','));

        // Background choice has to reach the chart panel and the paper.
        const doc = new DOMParser().parseFromString(svgStr, 'image/svg+xml');
        doc.querySelectorAll('.bglayer > rect.bg').forEach(r => {
            r.setAttribute('style', background === 'transparent' ? 'fill: none; stroke-width: 0;' : 'fill: #ffffff; fill-opacity: 1; stroke-width: 0;');
        });
        if (background === 'transparent') {
            const rw = doc.documentElement.getAttribute('width'), rh = doc.documentElement.getAttribute('height');
            doc.querySelectorAll('rect:not([class])').forEach(r => {
                if (r.getAttribute('width') === rw && r.getAttribute('height') === rh
                    && /fill:\s*rgb\(255,\s*255,\s*255\)/.test(r.getAttribute('style') || '')) r.setAttribute('style', 'fill: none;');
            });
        }
        svgStr = new XMLSerializer().serializeToString(doc.documentElement);
        await this._exportSvgString(svgStr, dlg, { filename, widthPx: w, heightPx: h });
    },

    // The dialog. Resolves { format, widthCm, heightCm, dpi, background } or null.
    // Default 10 cm wide, height from the on-screen aspect, 600 dpi; the last
    // choice is remembered for the session so a series of figures come out alike.
    _showExportDialog(context) {
        return new Promise(resolve => {
            const modal = document.getElementById('exportOptionsModal');
            if (!modal) { resolve(null); return; }
            const fmtEl = document.getElementById('exportOptFormat');
            const dpiRow = document.getElementById('exportOptDpiRow');
            const widthEl = document.getElementById('exportOptWidth');
            const heightEl = document.getElementById('exportOptHeight');
            const dpiEl = document.getElementById('exportOptDpi');
            const lockEl = document.getElementById('exportOptLockAspect');
            const aspect = context.plotH > 0 ? context.plotH / context.plotW : 1;
            this._exportDialogAspect = aspect;
            const prev = this._lastExportOpts || {};
            const defaultW = prev.widthCm || 10;
            widthEl.value = defaultW;
            heightEl.value = prev.lockAspect === false && prev.heightCm ? prev.heightCm : Math.max(2, Math.round(defaultW * aspect * 10) / 10);
            dpiEl.value = prev.dpi || 600;
            fmtEl.value = prev.format || 'png';
            lockEl.checked = prev.lockAspect !== false;
            document.querySelectorAll('input[name="exportOptBg"]').forEach(r => { r.checked = r.value === (prev.background || 'white'); });
            const syncFormat = () => { dpiRow.style.display = fmtEl.value === 'svg' || fmtEl.value === 'pdf' || fmtEl.value === 'pptx' ? 'none' : ''; };
            fmtEl.onchange = syncFormat; syncFormat();
            modal.style.display = 'flex';
            const cleanup = () => { modal.style.display = 'none'; };
            document.getElementById('exportOptConfirm').onclick = () => {
                const opts = {
                    format: fmtEl.value || 'png',
                    widthCm: this.numInput(widthEl, 10) || 10,
                    heightCm: this.numInput(heightEl, 10) || 10,
                    dpi: parseInt(dpiEl.value) || 600,
                    background: document.querySelector('input[name="exportOptBg"]:checked')?.value || 'white',
                    lockAspect: !!lockEl.checked
                };
                this._lastExportOpts = opts;
                cleanup(); resolve(opts);
            };
            const cancel = () => { cleanup(); resolve(null); };
            document.getElementById('exportOptCancel').onclick = cancel;
            document.getElementById('exportOptionsClose').onclick = cancel;
            modal.onclick = (e) => { if (e.target === modal) cancel(); };
        });
    },

    numInput(el, fallback) {
        const v = parseFloat(String(el?.value ?? '').replace(',', '.'));
        return isFinite(v) ? v : fallback;
    },

    adjustNumber(id, delta) {
        const el = document.getElementById(id);
        if (!el) return;
        const min = parseFloat(el.min), max = parseFloat(el.max);
        let v = this.numInput(el, 0) + delta;
        if (isFinite(min)) v = Math.max(min, v);
        if (isFinite(max)) v = Math.min(max, v);
        el.value = Math.round(v * 10) / 10;
    },

    adjustExportSize(id, delta) {
        this.adjustNumber(id, delta);
        const lockEl = document.getElementById('exportOptLockAspect');
        if (!lockEl || !lockEl.checked) return;
        const aspect = this._exportDialogAspect || 1;
        const widthEl = document.getElementById('exportOptWidth');
        const heightEl = document.getElementById('exportOptHeight');
        if (id === 'exportOptWidth') heightEl.value = Math.max(2, Math.round(this.numInput(widthEl, 10) * aspect * 10) / 10);
        else if (id === 'exportOptHeight') widthEl.value = aspect > 0 ? Math.max(2, Math.round((this.numInput(heightEl, 10) / aspect) * 10) / 10) : this.numInput(heightEl, 10);
    },

    applyExportSizePreset(widthCm) {
        const widthEl = document.getElementById('exportOptWidth');
        const heightEl = document.getElementById('exportOptHeight');
        const aspect = this._exportDialogAspect || 1;
        widthEl.value = widthCm;
        heightEl.value = Math.max(2, Math.round(widthCm * aspect * 10) / 10);
    },

    _saveBlob(blob, name) {
        const a = document.createElement('a');
        a.href = URL.createObjectURL(blob);
        a.download = name;
        document.body.appendChild(a); a.click(); document.body.removeChild(a);
        setTimeout(() => URL.revokeObjectURL(a.href), 2000);
    },

    // A finished SVG string goes out in whatever format the dialog chose.
    async _exportSvgString(svgIn, dlg, opts) {
        const { filename } = opts || {};
        const { widthCm, heightCm, dpi, background } = dlg;
        const fmt = dlg.format || 'png';
        const stem = `${filename}_${exportStamp()}`;
        let outSvg = svgIn;
        if (background === 'white' && !/id="enrichExportBg"/.test(outSvg)) {
            outSvg = outSvg.replace(/(<svg[^>]*>)/, `$1<rect id="enrichExportBg" x="0" y="0" width="100%" height="100%" fill="white"/>`);
        }
        if (fmt === 'svg') {
            // Absolute size on the outer tag so the file opens at the asked print size.
            const sized = outSvg.replace(/<svg\b[^>]*>/, (tag) => tag
                .replace(/\swidth="[^"]*"/, ` width="${widthCm}cm"`)
                .replace(/\sheight="[^"]*"/, ` height="${heightCm}cm"`));
            this._saveBlob(new Blob([sized], { type: 'image/svg+xml;charset=utf-8' }), `${stem}.svg`);
            return;
        }
        const targetPxW = Math.round(widthCm * dpi / 2.54);
        const targetPxH = Math.round(heightCm * dpi / 2.54);
        return new Promise(resolve => {
            const svgUrl = URL.createObjectURL(new Blob([outSvg], { type: 'image/svg+xml;charset=utf-8' }));
            const img = new Image();
            img.onload = async () => {
                const srcAR = (img.naturalHeight || targetPxH) / (img.naturalWidth || targetPxW);
                const drawW = targetPxW;
                const drawH = Math.abs(targetPxH / targetPxW - srcAR) > 0.01 ? Math.round(targetPxW * srcAR) : targetPxH;
                const canvas = document.createElement('canvas');
                canvas.width = drawW; canvas.height = drawH;
                const ctx = canvas.getContext('2d');
                if (background === 'white') { ctx.fillStyle = 'white'; ctx.fillRect(0, 0, drawW, drawH); }
                ctx.drawImage(img, 0, 0, drawW, drawH);
                URL.revokeObjectURL(svgUrl);
                const outH = Math.round(widthCm * (drawH / drawW) * 100) / 100;
                try {
                    await this._downloadCanvasAs(canvas, fmt, stem, { dpi, widthCm, heightCm: outH, svg: outSvg });
                } catch (e) { console.error('Export failed:', e); alert('Export failed: ' + e.message); }
                resolve();
            };
            img.onerror = () => { URL.revokeObjectURL(svgUrl); alert('Could not render the figure for export.'); resolve(); };
            img.src = svgUrl;
        });
    },

    async _downloadCanvasAs(canvas, fmt, stem, opts = {}) {
        const { dpi = 600, widthCm, heightCm, svg } = opts;
        const effDpi = (widthCm && canvas?.width) ? Math.max(1, Math.round(canvas.width / (widthCm / 2.54))) : dpi;
        if (fmt === 'tiff') { this._saveBlob(new Blob([this._canvasToTiff(canvas, effDpi)], { type: 'image/tiff' }), `${stem}.tiff`); return; }
        if (fmt === 'pdf') {
            if (svg && this._canVectorExport()) { this._saveBlob(await this._svgToPdfVector(svg, widthCm, heightCm), `${stem}.pdf`); return; }
            this._saveBlob(new Blob([this._canvasToPdf(canvas, widthCm, heightCm)], { type: 'application/pdf' }), `${stem}.pdf`); return;
        }
        if (fmt === 'pptx') { this._saveBlob(await this._canvasToPptx(canvas, widthCm, heightCm, svg), `${stem}.pptx`); return; }
        const durl = canvas.toDataURL('image/png');
        if (!durl || durl.length < 100) { alert('The image is too large at this size and DPI. Lower the DPI or the width.'); return; }
        let buf = await (await fetch(durl)).arrayBuffer();
        buf = this._setPngDpi(buf, effDpi);
        this._saveBlob(new Blob([buf], { type: 'image/png' }), `${stem}.png`);
    },

    // Baseline uncompressed RGB TIFF with the DPI in the resolution tags.
    _canvasToTiff(canvas, dpi) {
        const w = canvas.width, h = canvas.height;
        const data = canvas.getContext('2d').getImageData(0, 0, w, h).data;
        const strip = new Uint8Array(w * h * 3);
        for (let i = 0, j = 0; i < data.length; i += 4) {
            const a = data[i + 3] / 255;
            strip[j++] = Math.round(data[i] * a + 255 * (1 - a));
            strip[j++] = Math.round(data[i + 1] * a + 255 * (1 - a));
            strip[j++] = Math.round(data[i + 2] * a + 255 * (1 - a));
        }
        const nTags = 12, ifdSize = 2 + nTags * 12 + 4, extra = 8 + ifdSize;
        const bpsOff = extra, xresOff = extra + 6, yresOff = extra + 14, stripOff = extra + 22;
        const buf = new ArrayBuffer(stripOff + strip.length);
        const dv = new DataView(buf);
        dv.setUint16(0, 0x4949, true); dv.setUint16(2, 42, true); dv.setUint32(4, 8, true);
        let p = 8;
        dv.setUint16(p, nTags, true); p += 2;
        const tag = (id, type, count, value) => { dv.setUint16(p, id, true); dv.setUint16(p + 2, type, true); dv.setUint32(p + 4, count, true); dv.setUint32(p + 8, value, true); p += 12; };
        tag(256, 4, 1, w); tag(257, 4, 1, h); tag(258, 3, 3, bpsOff); tag(259, 3, 1, 1); tag(262, 3, 1, 2);
        tag(273, 4, 1, stripOff); tag(277, 3, 1, 3); tag(278, 4, 1, h); tag(279, 4, 1, strip.length);
        tag(282, 5, 1, xresOff); tag(283, 5, 1, yresOff); tag(296, 3, 1, 2);
        dv.setUint32(p, 0, true);
        dv.setUint16(bpsOff, 8, true); dv.setUint16(bpsOff + 2, 8, true); dv.setUint16(bpsOff + 4, 8, true);
        const d = Math.round(dpi) || 300;
        dv.setUint32(xresOff, d, true); dv.setUint32(xresOff + 4, 1, true);
        dv.setUint32(yresOff, d, true); dv.setUint32(yresOff + 4, 1, true);
        new Uint8Array(buf).set(strip, stripOff);
        return buf;
    },

    _canVectorExport() {
        const JS = window.jspdf?.jsPDF || window.jsPDF;
        return !!(JS && JS.API && typeof JS.API.svg === 'function');
    },

    // Vector PDF through jsPDF + svg2pdf, page sized to the asked cm.
    async _svgToPdfVector(svgStr, widthCm, heightCm) {
        const JS = window.jspdf?.jsPDF || window.jsPDF;
        if (!JS) throw new Error('jsPDF unavailable');
        // The standard PDF fonts lack the Unicode minus and most superscripts.
        svgStr = String(svgStr).replace(/−/g, '-');
        const SUP = { '⁻': '-', '⁰': '0', '¹': '1', '²': '2', '³': '3', '⁴': '4', '⁵': '5', '⁶': '6', '⁷': '7', '⁸': '8', '⁹': '9' };
        svgStr = svgStr.replace(/[⁻⁰¹²³⁴⁵⁶⁷⁸⁹]+/g, run => /[⁻⁰⁴⁵⁶⁷⁸⁹]/.test(run) ? '^' + [...run].map(c => SUP[c]).join('') : run);
        const ptW = (widthCm || 10) / 2.54 * 72, ptH = (heightCm || 10) / 2.54 * 72;
        const pdf = new JS({ unit: 'pt', orientation: ptW >= ptH ? 'landscape' : 'portrait', format: [Math.min(ptW, ptH), Math.max(ptW, ptH)], compress: true });
        const pw = pdf.internal.pageSize.getWidth(), ph = pdf.internal.pageSize.getHeight();
        const holder = document.createElement('div');
        holder.style.cssText = 'position:fixed; left:-99999px; top:0; width:0; height:0; overflow:hidden;';
        holder.innerHTML = svgStr;
        const svgEl = holder.querySelector('svg');
        if (!svgEl) throw new Error('no svg element');
        document.body.appendChild(holder);
        try { await pdf.svg(svgEl, { x: 0, y: 0, width: pw, height: ph }); }
        finally { document.body.removeChild(holder); }
        return pdf.output('blob');
    },

    // Raster fallback PDF: one page with the figure as JPEG.
    _canvasToPdf(canvas, widthCm, heightCm) {
        const tmp = document.createElement('canvas');
        tmp.width = canvas.width; tmp.height = canvas.height;
        const tctx = tmp.getContext('2d');
        tctx.fillStyle = '#fff'; tctx.fillRect(0, 0, tmp.width, tmp.height); tctx.drawImage(canvas, 0, 0);
        const bin = atob(tmp.toDataURL('image/jpeg', 0.95).split(',')[1]);
        const jpeg = new Uint8Array(bin.length);
        for (let i = 0; i < bin.length; i++) jpeg[i] = bin.charCodeAt(i);
        const ptW = (widthCm || 10) / 2.54 * 72, ptH = (heightCm || 10) / 2.54 * 72;
        const enc = new TextEncoder(); const parts = []; const offsets = []; let len = 0;
        const push = (s) => { const b = (typeof s === 'string') ? enc.encode(s) : s; parts.push(b); len += b.length; };
        push('%PDF-1.4\n');
        offsets.push(len); push('1 0 obj\n<< /Type /Catalog /Pages 2 0 R >>\nendobj\n');
        offsets.push(len); push('2 0 obj\n<< /Type /Pages /Kids [3 0 R] /Count 1 >>\nendobj\n');
        offsets.push(len); push(`3 0 obj\n<< /Type /Page /Parent 2 0 R /MediaBox [0 0 ${ptW.toFixed(2)} ${ptH.toFixed(2)}] /Resources << /XObject << /Im0 4 0 R >> >> /Contents 5 0 R >>\nendobj\n`);
        offsets.push(len); push(`4 0 obj\n<< /Type /XObject /Subtype /Image /Width ${canvas.width} /Height ${canvas.height} /ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /DCTDecode /Length ${jpeg.length} >>\nstream\n`);
        push(jpeg); push('\nendstream\nendobj\n');
        const content = `q ${ptW.toFixed(2)} 0 0 ${ptH.toFixed(2)} 0 0 cm /Im0 Do Q\n`;
        offsets.push(len); push(`5 0 obj\n<< /Length ${content.length} >>\nstream\n${content}endstream\nendobj\n`);
        const xrefStart = len;
        let xref = `xref\n0 6\n0000000000 65535 f \n`;
        for (const off of offsets) xref += String(off).padStart(10, '0') + ' 00000 n \n';
        push(xref); push(`trailer\n<< /Size 6 /Root 1 0 R >>\nstartxref\n${xrefStart}\n%%EOF`);
        const out = new Uint8Array(len); let o = 0;
        for (const pt of parts) { out.set(pt, o); o += pt.length; }
        return out.buffer;
    },

    // One 16:9 slide with the figure centred, as vector SVG with a PNG fallback.
    async _canvasToPptx(canvas, widthCm, heightCm, svgStr) {
        if (typeof JSZip === 'undefined') throw new Error('JSZip unavailable');
        const EMU = 360000, cx = 12192000, cy = 6858000;
        const figW = (widthCm || 10) * EMU, figH = (heightCm || 10) * EMU;
        const scale = Math.min(cx * 0.88 / figW, cy * 0.88 / figH);
        const picW = Math.round(figW * scale), picH = Math.round(figH * scale);
        const offX = Math.round((cx - picW) / 2), offY = Math.round((cy - picH) / 2);
        const pngB64 = canvas.toDataURL('image/png').split(',')[1];
        const useSvg = !!svgStr;
        const REL = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships';
        const zip = new JSZip();
        zip.file('[Content_Types].xml', `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types"><Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/><Default Extension="xml" ContentType="application/xml"/><Default Extension="png" ContentType="image/png"/>${useSvg ? '<Default Extension="svg" ContentType="image/svg+xml"/>' : ''}<Override PartName="/ppt/presentation.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.presentation.main+xml"/><Override PartName="/ppt/slideMasters/slideMaster1.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slideMaster+xml"/><Override PartName="/ppt/slideLayouts/slideLayout1.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slideLayout+xml"/><Override PartName="/ppt/slides/slide1.xml" ContentType="application/vnd.openxmlformats-officedocument.presentationml.slide+xml"/><Override PartName="/ppt/theme/theme1.xml" ContentType="application/vnd.openxmlformats-officedocument.theme+xml"/></Types>`);
        zip.file('_rels/.rels', `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="${REL}/officeDocument" Target="ppt/presentation.xml"/></Relationships>`);
        zip.file('ppt/presentation.xml', `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<p:presentation xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:r="${REL}" xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main"><p:sldMasterIdLst><p:sldMasterId id="2147483648" r:id="rId1"/></p:sldMasterIdLst><p:sldIdLst><p:sldId id="256" r:id="rId2"/></p:sldIdLst><p:sldSz cx="${cx}" cy="${cy}"/><p:notesSz cx="6858000" cy="9144000"/></p:presentation>`);
        zip.file('ppt/_rels/presentation.xml.rels', `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="${REL}/slideMaster" Target="slideMasters/slideMaster1.xml"/><Relationship Id="rId2" Type="${REL}/slide" Target="slides/slide1.xml"/><Relationship Id="rId3" Type="${REL}/theme" Target="theme/theme1.xml"/></Relationships>`);
        const clrMap = `<p:clrMap bg1="lt1" tx1="dk1" bg2="lt2" tx2="dk2" accent1="accent1" accent2="accent2" accent3="accent3" accent4="accent4" accent5="accent5" accent6="accent6" hlink="hlink" folHlink="folHlink"/>`;
        const emptyTree = `<p:spTree><p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr><p:grpSpPr/></p:spTree>`;
        zip.file('ppt/slideMasters/slideMaster1.xml', `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<p:sldMaster xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:r="${REL}" xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main"><p:cSld>${emptyTree}</p:cSld>${clrMap}<p:sldLayoutIdLst><p:sldLayoutId id="2147483649" r:id="rId1"/></p:sldLayoutIdLst></p:sldMaster>`);
        zip.file('ppt/slideMasters/_rels/slideMaster1.xml.rels', `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="${REL}/slideLayout" Target="../slideLayouts/slideLayout1.xml"/><Relationship Id="rId2" Type="${REL}/theme" Target="../theme/theme1.xml"/></Relationships>`);
        zip.file('ppt/slideLayouts/slideLayout1.xml', `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<p:sldLayout xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:r="${REL}" xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main" type="blank" preserve="1"><p:cSld name="Blank">${emptyTree}</p:cSld><p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:sldLayout>`);
        zip.file('ppt/slideLayouts/_rels/slideLayout1.xml.rels', `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="${REL}/slideMaster" Target="../slideMasters/slideMaster1.xml"/></Relationships>`);
        zip.file('ppt/theme/theme1.xml', this._minimalPptxTheme());
        const blip = useSvg
            ? `<a:blip r:embed="rId1"><a:extLst><a:ext uri="{96DAC541-7B7A-43D3-8B79-37D633B846F1}"><asvg:svgBlip xmlns:asvg="http://schemas.microsoft.com/office/drawing/2016/SVG/main" r:embed="rId3"/></a:ext></a:extLst></a:blip>`
            : `<a:blip r:embed="rId1"/>`;
        const pic = `<p:pic><p:nvPicPr><p:cNvPr id="2" name="Figure"/><p:cNvPicPr><a:picLocks noChangeAspect="1"/></p:cNvPicPr><p:nvPr/></p:nvPicPr><p:blipFill>${blip}<a:stretch><a:fillRect/></a:stretch></p:blipFill><p:spPr><a:xfrm><a:off x="${offX}" y="${offY}"/><a:ext cx="${picW}" cy="${picH}"/></a:xfrm><a:prstGeom prst="rect"><a:avLst/></a:prstGeom></p:spPr></p:pic>`;
        zip.file('ppt/slides/slide1.xml', `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<p:sld xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:r="${REL}" xmlns:p="http://schemas.openxmlformats.org/presentationml/2006/main"><p:cSld><p:spTree><p:nvGrpSpPr><p:cNvPr id="1" name=""/><p:cNvGrpSpPr/><p:nvPr/></p:nvGrpSpPr><p:grpSpPr/>${pic}</p:spTree></p:cSld><p:clrMapOvr><a:masterClrMapping/></p:clrMapOvr></p:sld>`);
        zip.file('ppt/slides/_rels/slide1.xml.rels', `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"><Relationship Id="rId1" Type="${REL}/image" Target="../media/image1.png"/><Relationship Id="rId2" Type="${REL}/slideLayout" Target="../slideLayouts/slideLayout1.xml"/>${useSvg ? `<Relationship Id="rId3" Type="${REL}/image" Target="../media/image2.svg"/>` : ''}</Relationships>`);
        zip.file('ppt/media/image1.png', pngB64, { base64: true });
        if (useSvg) zip.file('ppt/media/image2.svg', svgStr);
        return await zip.generateAsync({ type: 'blob', mimeType: 'application/vnd.openxmlformats-officedocument.presentationml.presentation' });
    },

    _minimalPptxTheme() {
        const A = 'http://schemas.openxmlformats.org/drawingml/2006/main';
        const solid = (c) => `<a:solidFill><a:srgbClr val="${c}"/></a:solidFill>`;
        const fillLst = `<a:fillStyleLst>${solid('FFFFFF')}${solid('FFFFFF')}${solid('FFFFFF')}</a:fillStyleLst>`;
        const lnLst = `<a:lnStyleLst><a:ln w="6350"><a:solidFill><a:srgbClr val="000000"/></a:solidFill></a:ln><a:ln w="12700"><a:solidFill><a:srgbClr val="000000"/></a:solidFill></a:ln><a:ln w="19050"><a:solidFill><a:srgbClr val="000000"/></a:solidFill></a:ln></a:lnStyleLst>`;
        const effLst = `<a:effectStyleLst><a:effectStyle><a:effectLst/></a:effectStyle><a:effectStyle><a:effectLst/></a:effectStyle><a:effectStyle><a:effectLst/></a:effectStyle></a:effectStyleLst>`;
        const bgLst = `<a:bgFillStyleLst>${solid('FFFFFF')}${solid('FFFFFF')}${solid('FFFFFF')}</a:bgFillStyleLst>`;
        const clr = `<a:clrScheme name="Office"><a:dk1><a:sysClr val="windowText" lastClr="000000"/></a:dk1><a:lt1><a:sysClr val="window" lastClr="FFFFFF"/></a:lt1><a:dk2><a:srgbClr val="44546A"/></a:dk2><a:lt2><a:srgbClr val="E7E6E6"/></a:lt2><a:accent1><a:srgbClr val="4472C4"/></a:accent1><a:accent2><a:srgbClr val="ED7D31"/></a:accent2><a:accent3><a:srgbClr val="A5A5A5"/></a:accent3><a:accent4><a:srgbClr val="FFC000"/></a:accent4><a:accent5><a:srgbClr val="5B9BD5"/></a:accent5><a:accent6><a:srgbClr val="70AD47"/></a:accent6><a:hlink><a:srgbClr val="0563C1"/></a:hlink><a:folHlink><a:srgbClr val="954F72"/></a:folHlink></a:clrScheme>`;
        const font = `<a:fontScheme name="Office"><a:majorFont><a:latin typeface="Calibri Light"/><a:ea typeface=""/><a:cs typeface=""/></a:majorFont><a:minorFont><a:latin typeface="Calibri"/><a:ea typeface=""/><a:cs typeface=""/></a:minorFont></a:fontScheme>`;
        return `<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n<a:theme xmlns:a="${A}" name="Office"><a:themeElements>${clr}${font}<a:fmtScheme name="Office">${fillLst}${lnLst}${effLst}${bgLst}</a:fmtScheme></a:themeElements></a:theme>`;
    },

    // Write a pHYs chunk so the PNG carries its DPI.
    _setPngDpi(arrayBuffer, dpi) {
        const ppm = Math.round(dpi * 39.3701);
        const src = new Uint8Array(arrayBuffer);
        const T = new Uint32Array(256);
        for (let n = 0; n < 256; n++) { let c = n; for (let k = 0; k < 8; k++) c = (c & 1) ? (0xedb88320 ^ (c >>> 1)) : (c >>> 1); T[n] = c >>> 0; }
        const crc32 = (bytes) => { let c = 0xffffffff; for (let i = 0; i < bytes.length; i++) c = T[(c ^ bytes[i]) & 0xff] ^ (c >>> 8); return (c ^ 0xffffffff) >>> 0; };
        const physData = new Uint8Array(9);
        const dv = new DataView(physData.buffer);
        dv.setUint32(0, ppm); dv.setUint32(4, ppm); physData[8] = 1;
        const typeAndData = new Uint8Array(13);
        typeAndData.set([0x70, 0x48, 0x59, 0x73]); typeAndData.set(physData, 4);
        const phys = new Uint8Array(21);
        new DataView(phys.buffer).setUint32(0, 9); phys.set(typeAndData, 4);
        new DataView(phys.buffer).setUint32(17, crc32(typeAndData));
        let pos = 8; const chunks = [];
        while (pos < src.length) {
            const len = new DataView(src.buffer, src.byteOffset + pos, 4).getUint32(0);
            const type = String.fromCharCode(src[pos + 4], src[pos + 5], src[pos + 6], src[pos + 7]);
            const chunkLen = 12 + len;
            if (type !== 'pHYs') chunks.push(src.subarray(pos, pos + chunkLen));
            pos += chunkLen;
        }
        const out = [src.subarray(0, 8)]; let inserted = false;
        for (const c of chunks) { out.push(c); if (!inserted && String.fromCharCode(c[4], c[5], c[6], c[7]) === 'IHDR') { out.push(phys); inserted = true; } }
        let total = 0; for (const c of out) total += c.length;
        const result = new Uint8Array(total); let off = 0;
        for (const c of out) { result.set(c, off); off += c.length; }
        return result.buffer;
    }
});
