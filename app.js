/* ============================================================
   GSE332551 portal app logic
   ============================================================ */

/* ---------- provenance: commit stamp + per-download SHA-256 (A4/A5) ---------- */
(() => {
  const b = window.NF_BUILD || {};
  const c = document.getElementById("nf-commit");
  const d = document.getElementById("nf-date");
  if (c) {
    const sha = b.sha && b.sha !== "dev" ? b.sha : "local";
    c.innerHTML = `<code>${sha}</code>`;
    if (b.commit) { c.href = b.commit; c.target = "_blank"; c.rel = "noopener"; }
  }
  if (d) d.textContent = b.date || "unbuilt";
})();

// Annotate every download link with the SHA-256 of the file it points to.
fetch("./downloads/checksums.json")
  .then((r) => (r.ok ? r.json() : {}))
  .then((sums) => {
    const base = (href) => (href || "").split("/").pop().split("?")[0];
    document.querySelectorAll("a[download]").forEach((a) => {
      const sha = sums[base(a.getAttribute("href"))];
      if (!sha) return;
      a.title = `sha256 ${sha}`;
      if (a.classList.contains("dl-bulk-card")) {
        const s = document.createElement("div");
        s.className = "dbc-sha";
        s.textContent = `sha256 ${sha.slice(0, 20)}…`;
        a.appendChild(s);
      } else if (a.closest(".dl-row")) {
        const s = document.createElement("div");
        s.className = "dl-sha";
        s.textContent = `sha256 ${sha}`;
        a.closest(".dl-row").querySelector(".info").appendChild(s);
      }
    });
    // Figure captions: show the SHA-256 of the displayed figure file.
    document.querySelectorAll("figure.fig").forEach((fig) => {
      const img = fig.querySelector("img");
      const cap = fig.querySelector("figcaption");
      const sha = img && sums[base(img.getAttribute("src"))];
      if (!sha || !cap) return;
      const s = document.createElement("div");
      s.className = "fig-prov";
      s.textContent = `sha256 ${sha}`;
      cap.appendChild(s);
    });
  })
  .catch(() => {});

/* ---------- theme ---------- */
const root = document.documentElement;
const savedTheme = localStorage.getItem("nf-theme");
if (savedTheme) root.setAttribute("data-theme", savedTheme);
document.getElementById("theme-toggle").addEventListener("click", () => {
  const next = root.getAttribute("data-theme") === "dark" ? "light" : "dark";
  root.setAttribute("data-theme", next);
  localStorage.setItem("nf-theme", next);
});

/* ---------- mobile nav ---------- */
const toc = document.getElementById("toc");
document.getElementById("hamburger").addEventListener("click", () => toc.classList.toggle("open"));
toc.querySelectorAll("a").forEach((a) => a.addEventListener("click", () => toc.classList.remove("open")));

/* ---------- lock (clear the gate cookie, reload) ---------- */
document.getElementById("lock").addEventListener("click", () => {
  document.cookie = "nf_gate=; Path=/; Max-Age=0; SameSite=Lax";
  location.reload();
});

/* ---------- scroll-spy TOC ---------- */
const tocLinks = [...toc.querySelectorAll("a")];
const spy = new IntersectionObserver(
  (entries) => {
    entries.forEach((e) => {
      if (e.isIntersecting) {
        const id = e.target.id;
        tocLinks.forEach((l) => l.classList.toggle("active", l.getAttribute("href") === `#${id}`));
      }
    });
  },
  { rootMargin: "-40% 0px -55% 0px" }
);
document.querySelectorAll("main .section").forEach((s) => spy.observe(s));

/* ---------- lightbox ---------- */
const lb = document.getElementById("lightbox");
const lbImg = document.getElementById("lightbox-img");
document.querySelectorAll(".fig img").forEach((img) => {
  img.addEventListener("click", () => {
    lbImg.src = img.src;
    lb.hidden = false;
  });
});
function closeLb() { lb.hidden = true; lbImg.src = ""; }
document.getElementById("lightbox-close").addEventListener("click", closeLb);
lb.addEventListener("click", (e) => { if (e.target === lb) closeLb(); });
document.addEventListener("keydown", (e) => { if (e.key === "Escape") closeLb(); });

/* ---------- scripts list ---------- */
const SCRIPTS = [
  ["01_deseq2_rnaseq.R", "r", "Original RNA-seq DESeq2: transcripts collapsed to genes, filter mean ≥ 10, AD vs CN."],
  ["02_deseq2_sncrnaseq.R", "r", "Small-RNA DESeq2 run directly on sncRNA features; classifies piRNA, tRF, miRNA, snoRNA."],
  ["03_visualize.py", "py", "Volcano plots, sncRNA class breakdown, and heatmaps for the original run."],
  ["04_paper_comparison.py", "py", "Gene-by-gene and sncRNA-by-sncRNA comparison of reanalysis vs the paper."],
  ["05_reproduce_paper.R", "r", "Transcript-level (paper approach) and gene-level side by side to diagnose the 827-vs-195 gap."],
  ["06_diagnose_normalization.R", "r", "Library sizes, size factors, correlations, PCA. Identifies R17_AD as the outlier."],
  ["07_corrected_DEGs.R", "r", "Definitive corrected DEG list: gene-level DESeq2, R17_AD excluded, LFC shrinkage, FDR."],
  ["08_visualize_corrected.py", "py", "Volcano and high-confidence DEG heatmap from the corrected analysis."],
  ["09_gsea_pathways.py", "py", "Preranked GSEA (Wald-ranked) vs Hallmark, GO-BP, Reactome, KEGG."],
  ["10_sncrna_de.R", "r", "Small-RNA DESeq2, classified into miRNA / piRNA / tRF / snoRNA."],
  ["11_mirna_mrna_integration.py", "py", "Cross-assay test: do DE-miRNA targets move oppositely in the RNA-seq?"],
  ["12_microglia_states.py", "py", "ssGSEA scoring of microglial states (homeostatic / DAM / LDAM / SREBP / NF-kB)."],
  ["13_tf_activity.py", "py", "TF-activity inference (decoupler ULM + CollecTRI) on the AD-vs-CN contrast."],
];
document.getElementById("script-list").innerHTML = SCRIPTS.map(([name, lang, desc], i) => `
  <div class="dl-row">
    <div class="dl-ico">${String(i + 1).padStart(2, "0")}</div>
    <div class="info">
      <div class="fname">${name}<span class="lang ${lang}">${lang === "r" ? "R" : "Python"}</span></div>
      <div class="desc">${desc}</div>
    </div>
    <button class="btn btn-sm ide-open" data-i="${i}">&lt;/&gt; View code</button>
  </div>`).join("");

/* ---------- VSCode-style code viewer (Monaco Editor via CDN, no backend) ---------- */
const MONACO_BASE = "https://cdn.jsdelivr.net/npm/monaco-editor@0.52.2/min/vs";
const IDE_FILES = SCRIPTS.map(([name, lang]) => ({ name, lang, model: null }));
let monacoReady = null;
let ideEditor = null;

function loadMonaco() {
  if (monacoReady) return monacoReady;
  // route Monaco's workers through a same-context proxy so the CDN load works
  window.MonacoEnvironment = {
    getWorkerUrl: () =>
      "data:text/javascript;charset=utf-8," +
      encodeURIComponent(
        "self.MonacoEnvironment={baseUrl:'https://cdn.jsdelivr.net/npm/monaco-editor@0.52.2/min/'};" +
        `importScripts('${MONACO_BASE}/base/worker/workerMain.js');`
      ),
  };
  monacoReady = new Promise((resolve, reject) => {
    const s = document.createElement("script");
    s.src = `${MONACO_BASE}/loader.js`;
    s.onload = () => {
      window.require.config({ paths: { vs: MONACO_BASE } });
      window.require(["vs/editor/editor.main"], () => resolve(window.monaco));
    };
    s.onerror = () => reject(new Error("Monaco failed to load"));
    document.head.appendChild(s);
  });
  return monacoReady;
}

const ideOverlay = document.getElementById("code-ide");
const ideFilesEl = document.getElementById("ide-files");
const ideTab = document.getElementById("ide-tab");
const ideLang = document.getElementById("ide-lang");
const badge = (l) => `<span class="fx ${l}">${l === "r" ? "R" : "PY"}</span>`;

ideFilesEl.innerHTML = IDE_FILES.map((f, i) =>
  `<div class="ide-file" data-i="${i}">${badge(f.lang)}${f.name}</div>`).join("");

async function showFile(i) {
  const monaco = await loadMonaco();
  const f = IDE_FILES[i];
  if (!f.model) {
    let text;
    try {
      text = await (await fetch(`./scripts/${f.name}`)).text();
    } catch {
      text = `# Could not load ${f.name}`;
    }
    f.model = monaco.editor.createModel(text, f.lang === "r" ? "r" : "python");
  }
  ideEditor.setModel(f.model);
  ideEditor.setScrollPosition({ scrollTop: 0 });
  ideTab.innerHTML = `${badge(f.lang)} ${f.name}`;
  ideLang.textContent = f.lang === "r" ? "R" : "Python";
  [...ideFilesEl.children].forEach((el, idx) => el.classList.toggle("active", idx === i));
}

async function openIDE(i = 0) {
  ideOverlay.hidden = false;
  document.body.style.overflow = "hidden";
  let monaco;
  try {
    monaco = await loadMonaco();
  } catch {
    ideTab.textContent = "Editor failed to load (network offline?)";
    return;
  }
  if (!ideEditor) {
    ideEditor = monaco.editor.create(document.getElementById("monaco"), {
      value: "",
      language: "python",
      theme: "vs-dark",
      readOnly: true,
      automaticLayout: true,
      fontSize: 13,
      minimap: { enabled: true },
      scrollBeyondLastLine: false,
      renderWhitespace: "none",
    });
  }
  showFile(i);
}

function closeIDE() {
  ideOverlay.hidden = true;
  document.body.style.overflow = "";
}

document.getElementById("open-code").addEventListener("click", () => openIDE(0));
document.getElementById("ide-close").addEventListener("click", closeIDE);
ideFilesEl.addEventListener("click", (e) => {
  const row = e.target.closest(".ide-file");
  if (row) showFile(Number(row.dataset.i));
});
document.getElementById("script-list").addEventListener("click", (e) => {
  const btn = e.target.closest(".ide-open");
  if (btn) openIDE(Number(btn.dataset.i));
});
document.addEventListener("keydown", (e) => {
  if (e.key === "Escape" && !ideOverlay.hidden) closeIDE();
});

/* ---------- helpers ---------- */
// Link a gene symbol to its locus in the UCSC Genome Browser (GRCh38/hg38).
const GENOME_DB = "hg38";
function geneCell(sym) {
  const url = `https://genome.ucsc.edu/cgi-bin/hgTracks?db=${GENOME_DB}&position=${encodeURIComponent(sym)}`;
  return `<a class="gene-link" href="${url}" target="_blank" rel="noopener"` +
    ` title="Open ${sym} at its locus in the UCSC Genome Browser (${GENOME_DB})"><strong>${sym}</strong></a>`;
}

function fmtP(v) {
  if (v === null || v === undefined) return "n/a";
  if (v === 0) return "0";
  return v < 1e-3 ? v.toExponential(2) : v.toFixed(4);
}
function makeSortable(theadRow, getState, render) {
  theadRow.querySelectorAll("th[data-sort]").forEach((th) => {
    th.addEventListener("click", () => {
      const s = getState();
      const key = th.dataset.sort;
      if (s.key === key) s.dir *= -1;
      else { s.key = key; s.dir = key === "gene" || key === "pvalue" || key === "padj" || key === "fdr" ? 1 : -1; }
      render();
    });
  });
}

/* ---------- high-confidence genes table ---------- */
const GENES = [
  { gene: "PHLDA1", raw: 2.059, shrunk: 1.495, fdr: 4.51e-11, dir: "up_in_AD" },
  { gene: "SLC44A2", raw: -1.363, shrunk: -0.988, fdr: 6.59e-4, dir: "down_in_AD" },
  { gene: "ALK", raw: 2.794, shrunk: 1.053, fdr: 1.32e-3, dir: "up_in_AD" },
  { gene: "KCNQ3", raw: -1.156, shrunk: -0.885, fdr: 1.77e-3, dir: "down_in_AD" },
  { gene: "LRRC4", raw: -1.318, shrunk: -0.927, fdr: 2.57e-3, dir: "down_in_AD" },
  { gene: "PLOD2", raw: -2.552, shrunk: -0.840, fdr: 8.07e-3, dir: "down_in_AD" },
  { gene: "CKB", raw: -2.185, shrunk: -0.873, fdr: 8.31e-3, dir: "down_in_AD" },
  { gene: "ZFP3", raw: -2.292, shrunk: -0.865, fdr: 8.31e-3, dir: "down_in_AD" },
  { gene: "LRP8", raw: 1.179, shrunk: 0.837, fdr: 1.52e-2, dir: "up_in_AD" },
  { gene: "RYR1", raw: -1.712, shrunk: -0.869, fdr: 1.77e-2, dir: "down_in_AD" },
  { gene: "AP1S2", raw: 1.406, shrunk: 0.868, fdr: 2.45e-2, dir: "up_in_AD" },
  { gene: "MID1IP1", raw: 1.270, shrunk: 0.825, fdr: 3.24e-2, dir: "up_in_AD" },
  { gene: "SLIT3", raw: 1.378, shrunk: 0.837, fdr: 3.44e-2, dir: "up_in_AD" },
];
const geneState = { key: "fdr", dir: 1 };
const geneTable = document.getElementById("gene-table");
const geneFilter = document.getElementById("gene-filter");
function renderGenes() {
  const q = geneFilter.value.trim().toLowerCase();
  let rows = GENES.filter((r) => `${r.gene} ${r.dir}`.toLowerCase().includes(q));
  const rank = (r) => (r.dir === "up_in_AD" ? 1 : -1);
  rows = rows.slice().sort((a, b) => {
    let av = geneState.key === "dir" ? rank(a) : a[geneState.key];
    let bv = geneState.key === "dir" ? rank(b) : b[geneState.key];
    if (typeof av === "string") return av < bv ? -geneState.dir : av > bv ? geneState.dir : 0;
    return (av - bv) * geneState.dir;
  });
  geneTable.innerHTML = rows.map((r) => {
    const cls = r.dir === "up_in_AD" ? "cell-up" : "cell-down";
    return `<tr>
      <td>${geneCell(r.gene)}</td>
      <td class="num">${r.raw.toFixed(3)}</td>
      <td class="num">${r.shrunk.toFixed(3)}</td>
      <td class="num">${fmtP(r.fdr)}</td>
      <td class="${cls}">${r.dir.replaceAll("_", " ")}</td>
    </tr>`;
  }).join("");
}
makeSortable(document.querySelector("#gene-table-wrap thead tr"), () => geneState, renderGenes);
geneFilter.addEventListener("input", renderGenes);
renderGenes();

/* ---------- interactive DEG explorer ---------- */
const DEG = (window.DEG_ALL || []).map(([gene, baseMean, raw, shrunk, pvalue, padj]) => ({
  gene, baseMean, raw, shrunk, pvalue, padj,
  dir: raw >= 0 ? "up_in_AD" : "down_in_AD",
}));
const pRange = document.getElementById("p-range");
const fdrRange = document.getElementById("fdr-range");
const fcRange = document.getElementById("fc-range");
const degSearch = document.getElementById("deg-search");
const pVal = document.getElementById("p-val");
const fdrVal = document.getElementById("fdr-val");
const fcVal = document.getElementById("fc-val");
const summary = document.getElementById("deg-summary");
const degBody = document.querySelector("#deg-explorer-table tbody");
const foot = document.getElementById("deg-foot");
const MAX_SHOWN = 400;
const degState = { key: "pvalue", dir: 1 };

function passes(row) {
  const pMax = parseFloat(pRange.value);
  const fdrMax = parseFloat(fdrRange.value);
  const fcMin = parseFloat(fcRange.value);
  if (row.pvalue > pMax) return false;
  if (fdrMax < 1 && (row.padj === null || row.padj > fdrMax)) return false;
  if (Math.abs(row.raw) < fcMin) return false;
  return true;
}
function renderExplorer() {
  if (!DEG.length) return;
  pVal.textContent = parseFloat(pRange.value).toFixed(3);
  fdrVal.textContent = parseFloat(fdrRange.value).toFixed(3);
  fcVal.textContent = parseFloat(fcRange.value).toFixed(2);

  const term = degSearch.value.trim().toLowerCase();
  const hits = DEG.filter(passes);
  const up = hits.filter((r) => r.dir === "up_in_AD").length;
  summary.innerHTML = `
    <span class="pill total">${hits.length} DEGs</span>
    <span class="pill up">${up} up in AD</span>
    <span class="pill down">${hits.length - up} down in AD</span>`;

  let shown = term ? hits.filter((r) => r.gene.toLowerCase().includes(term)) : hits;
  const rank = (r) => (r.dir === "up_in_AD" ? 1 : -1);
  shown = shown.slice().sort((a, b) => {
    let av, bv;
    if (degState.key === "gene") { av = a.gene; bv = b.gene; }
    else if (degState.key === "dir") { av = rank(a); bv = rank(b); }
    else { av = a[degState.key]; bv = b[degState.key]; if (av === null) av = Infinity; if (bv === null) bv = Infinity; }
    if (av < bv) return -degState.dir;
    if (av > bv) return degState.dir;
    return 0;
  });

  const clipped = shown.length > MAX_SHOWN;
  degBody.innerHTML = (clipped ? shown.slice(0, MAX_SHOWN) : shown).map((r) => {
    const cls = r.dir === "up_in_AD" ? "cell-up" : "cell-down";
    return `<tr>
      <td>${geneCell(r.gene)}</td>
      <td class="num">${r.baseMean.toLocaleString()}</td>
      <td class="num">${r.raw.toFixed(3)}</td>
      <td class="num">${r.shrunk.toFixed(3)}</td>
      <td class="num">${fmtP(r.pvalue)}</td>
      <td class="num">${fmtP(r.padj)}</td>
      <td class="${cls}">${r.dir.replaceAll("_", " ")}</td>
    </tr>`;
  }).join("");

  foot.textContent = clipped
    ? `Showing first ${MAX_SHOWN} of ${shown.length} matching genes (sorted). Narrow the thresholds or search to see the rest.`
    : `Showing ${shown.length} matching gene${shown.length === 1 ? "" : "s"}.`;
}
if (DEG.length) {
  [pRange, fdrRange, fcRange].forEach((el) => el.addEventListener("input", renderExplorer));
  degSearch.addEventListener("input", renderExplorer);
  makeSortable(document.querySelector("#deg-explorer-table thead tr"), () => degState, renderExplorer);
  renderExplorer();
}
