// Builds presentation/exlink_onepager.pptx — a one-slide overview of the study,
// with entrance animations injected into the slide XML after pptxgenjs writes it.
// Usage: node presentation/build_onepager.js   (from the repository root)
const path = require("path");
const fs = require("fs");
const pptxgen = require("pptxgenjs");
const React = require("react");
const ReactDOMServer = require("react-dom/server");
const sharp = require("sharp");
const JSZip = require("jszip");
const { FaGasPump, FaProjectDiagram, FaShieldAlt, FaExclamationTriangle, FaBalanceScale, FaCrosshairs } = require("react-icons/fa");

const ROOT = path.resolve(__dirname, "..");
const OUT = path.join(__dirname, "exlink_onepager.pptx");
const APPLY_THEME = process.env.APPLY_THEME_JS; // optional: path to the pptx skill's apply_theme.js

const THEME = {
  name: "EX-link",
  headFontFace: "Cambria",
  bodyFontFace: "Calibri",
  colors: {
    dk1: "14181D", lt1: "FFFFFF", dk2: "212932", lt2: "AAB6C3",
    accent1: "F4A259", accent2: "4FB3A9", accent3: "E76F51", accent4: "8AB17D",
    accent5: "6C8EAD", accent6: "C9CED6", hlink: "F4A259", folHlink: "4FB3A9",
  },
};

async function icon(Comp, color) {
  const svg = ReactDOMServer.renderToStaticMarkup(React.createElement(Comp, { color: "#" + color, size: 256 }));
  const buf = await sharp(Buffer.from(svg)).resize(256, 256).png().toBuffer();
  return "image/png;base64," + buf.toString("base64");
}

async function main() {
  const pres = new pptxgen();
  pres.layout = "LAYOUT_WIDE"; // 13.33" x 7.5"
  pres.title = "EX-link Atkinson-cycle engine optimization";
  pres.author = "Simone Coniglio";
  pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
  const C = pres.SchemeColor;

  pres.defineSlideMaster({
    title: "ONE_PAGER_DARK",
    background: { color: C.text1 },
    objects: [
      { placeholder: { options: { name: "title", type: "title", x: 0.5, y: 0.72, w: 6.2, h: 1.25,
        fontFace: THEME.headFontFace, fontSize: 36, bold: true, color: C.background1, valign: "top", align: "left", margin: 0 },
        text: "Title" } },
      { text: { text: "github.com/SimoneConiglio/Atkinson-cycle-engine-optimisation-",
        options: { x: 6.95, y: 7.0, w: 5.9, h: 0.3, fontSize: 10, color: C.background2, align: "right", margin: 0 } } },
    ],
  });

  pres.addSection({ title: "Overview" });
  const s = pres.addSlide({ masterName: "ONE_PAGER_DARK", sectionTitle: "Overview" });
  s.transition = undefined;

  // ---- header ----
  s.addText("SHELL ECO-MARATHON  ·  MULTIDISCIPLINARY DESIGN  ·  GEMSEO", {
    objectName: "kicker", isTextBox: true, x: 0.5, y: 0.4, w: 6.2, h: 0.3, margin: 0,
    fontSize: 12, bold: true, charSpacing: 2, color: C.accent1 });
  s.addText("Range-optimal Atkinson engine linkage", { placeholder: "title" });
  s.addText("An extended-expansion (EX-link) engine designed for range — distance per litre of fuel — pricing efficiency, envelope, torque ripple and structural mass at the rates the physics fixes.", {
    objectName: "subtitle", isTextBox: true, x: 0.5, y: 1.95, w: 6.2, h: 0.8, margin: 0,
    fontSize: 14, color: C.background2, valign: "top" });

  // ---- animated mechanism ----
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { objectName: "gif-frame", x: 6.85, y: 0.4, w: 6.0, h: 3.85,
    rectRadius: 0.12, fill: { color: C.background1 }, line: { color: C.background1, width: 0 } });
  s.addImage({ objectName: "mechanism-gif", path: path.join(ROOT, "docs/figures/exlink.gif"),
    x: 6.92, y: 0.455, w: 5.86, h: 3.6625,
    altText: "The optimised EX-link linkage turning through one cycle, with piston motion, pressure-volume cycle and crankshaft torque" });
  s.addText("The design the study arrives at, through one cycle: 360° of the half-speed shaft, 720° of the crankshaft.", {
    objectName: "gif-caption", isTextBox: true, x: 6.85, y: 4.33, w: 6.0, h: 0.4, margin: 0,
    fontSize: 10, italic: true, color: C.background2 });

  // ---- stat tiles (2 x 2) ----
  const stats = [
    { v: "3395", u: "km/L", l: "range at failure probability 1.3 × 10⁻³", c: C.accent1 },
    { v: "+17.6", u: "%", l: "over a conventional engine (2888 km/L)", c: C.accent2 },
    { v: "+19.3", u: "%", l: "over a four-speed schedule (3260 vs 2733)", c: C.accent2 },
    { v: "12.9", u: "kg", l: "engine mass, against 16.9 kg conventional", c: C.accent1 },
  ];
  const tileW = 3.0, tileH = 0.85, gap = 0.2, tx = 0.5, ty = 2.95;
  stats.forEach((st, i) => {
    const x = tx + (i % 2) * (tileW + gap), y = ty + Math.floor(i / 2) * (tileH + gap);
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { objectName: `stat-${i}-bg`, x, y, w: tileW, h: tileH,
      rectRadius: 0.08, fill: { color: C.text2 }, line: { color: C.text2, width: 0 } });
    s.addText([
      { text: st.v, options: { fontSize: 30, bold: true, color: st.c, fontFace: THEME.headFontFace } },
      { text: " " + st.u, options: { fontSize: 14, bold: true, color: st.c } },
    ], { objectName: `stat-${i}-value`, isTextBox: true, x: x + 0.18, y: y + 0.05, w: tileW - 0.3, h: 0.45, margin: 0, valign: "middle" });
    s.addText(st.l, { objectName: `stat-${i}-label`, isTextBox: true, x: x + 0.18, y: y + 0.5, w: tileW - 0.3, h: 0.28,
      margin: 0, fontSize: 11, color: C.background2, valign: "top" });
  });

  // ---- findings ----
  const iconColor = THEME.colors.dk1;
  const findings = [
    { ic: FaExclamationTriangle, c: C.accent3, h: "The quasi-static optimum is the worst place to be",
      b: "Ignoring inertia drives the linkage to its transmission-angle singularity: no feasible structure above 2000 rpm. Backed off, it weighs half as much and goes further." },
    { ic: FaBalanceScale, c: C.accent1, h: "The topology is worth 17.6 %",
      b: "Indicated efficiency 0.457 → 0.480, mechanical 0.787 → 0.865, mass 16.9 → 12.9 kg. Extended expansion is the smallest of the three." },
    { ic: FaShieldAlt, c: C.accent2, h: "Reliability decides which bounds are real",
      b: "The reference design misses a requirement 66.4 % of the time. Holding constraints and a reliability target during the search gives 3395 km/L at Pf = 10⁻³." },
  ];
  const fy = 5.45, fh = 1.45, fw = (12.33 - 2 * 0.3) / 3;
  for (let i = 0; i < findings.length; i++) {
    const f = findings[i];
    const x = 0.5 + i * (fw + 0.3);
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { objectName: `finding-${i}-bg`, x, y: fy, w: fw, h: fh,
      rectRadius: 0.08, fill: { color: C.text2 }, line: { color: C.text2, width: 0 } });
    s.addShape(pres.shapes.OVAL, { objectName: `finding-${i}-circle`, x: x + 0.2, y: fy + 0.2, w: 0.5, h: 0.5,
      fill: { color: f.c }, line: { color: f.c, width: 0 } });
    s.addImage({ objectName: `finding-${i}-icon`, data: await icon(f.ic, iconColor), x: x + 0.32, y: fy + 0.32, w: 0.26, h: 0.26,
      altText: f.h });
    s.addText(f.h, { objectName: `finding-${i}-head`, isTextBox: true, x: x + 0.85, y: fy + 0.17, w: fw - 1.0, h: 0.56,
      margin: 0, fontSize: 14, bold: true, color: C.background1, valign: "middle" });
    s.addText(f.b, { objectName: `finding-${i}-body`, isTextBox: true, x: x + 0.2, y: fy + 0.78, w: fw - 0.4, h: 0.6,
      margin: 0, fontSize: 11, color: C.background2, valign: "top" });
  }
  s.addText("Five coupled disciplines solved as MDF with a Gauss–Seidel MDA · exact analytic derivatives · FORM reliability over 17 uncertain parameters", {
    objectName: "method-line", isTextBox: true, x: 0.5, y: 5.03, w: 12.33, h: 0.3, margin: 0,
    fontSize: 12, italic: true, color: C.accent6 });

  s.addNotes("One-page overview of the EX-link Atkinson-cycle engine study. The headline: optimising for range, rather than the conventional lever-arm quality measure, gives 3395 km/L at a failure probability of 1.3e-3, 17.6 % further than a conventional engine sized by the same models. Three findings: the quasi-static optimum sits at a singularity and cannot be built above 2000 rpm; the topology advantage comes mostly from mechanical efficiency and flywheel mass, not extended expansion; and reliability analysis shows which stated bounds actually bind.");

  await pres.writeFile({ fileName: OUT });
  if (APPLY_THEME) await require(APPLY_THEME).applyTheme(OUT, THEME);
  await addAnimations(OUT);
  console.log("wrote", OUT);
}

// ---------------------------------------------------------------------------
// Entrance animations. pptxgenjs cannot write <p:timing>, so we add it here.
// Everything plays automatically on slide start, one group after another.
// ---------------------------------------------------------------------------
const SEQUENCE = [
  // [start ms, effect, duration ms, object names]
  [0,    "fade", 600, ["kicker"]],
  [300,  "fade", 700, ["title"]],
  [800,  "fade", 700, ["subtitle"]],
  [1200, "zoom", 800, ["gif-frame", "mechanism-gif"]],
  [1800, "fade", 500, ["gif-caption"]],
  [2000, "fly",  500, ["stat-0-bg", "stat-0-value", "stat-0-label"]],
  [2250, "fly",  500, ["stat-1-bg", "stat-1-value", "stat-1-label"]],
  [2500, "fly",  500, ["stat-2-bg", "stat-2-value", "stat-2-label"]],
  [2750, "fly",  500, ["stat-3-bg", "stat-3-value", "stat-3-label"]],
  [3300, "wipe", 600, ["method-line"]],
  [3800, "fly",  600, ["finding-0-bg", "finding-0-circle", "finding-0-icon", "finding-0-head", "finding-0-body"]],
  [4200, "fly",  600, ["finding-1-bg", "finding-1-circle", "finding-1-icon", "finding-1-head", "finding-1-body"]],
  [4600, "fly",  600, ["finding-2-bg", "finding-2-circle", "finding-2-icon", "finding-2-head", "finding-2-body"]],
];

async function addAnimations(file) {
  const zip = await JSZip.loadAsync(fs.readFileSync(file));
  const slidePath = "ppt/slides/slide1.xml";
  let xml = await zip.file(slidePath).async("string");

  const ids = {};
  for (const m of xml.matchAll(/<p:cNvPr id="(\d+)" name="([^"]+)"/g)) ids[m[2]] = m[1];
  // the title placeholder is named by pptxgenjs, not by objectName
  const titleMatch = xml.match(/<p:cNvPr id="(\d+)" name="[^"]*"[^>]*>(?:(?!<p:cNvPr).)*?<p:ph[^>]*type="title"/s);
  if (titleMatch) ids.title = titleMatch[1];

  let n = 3;
  const tgt = (id) => `<p:tgtEl><p:spTgt spid="${id}"/></p:tgtEl>`;
  const setVisible = (id) =>
    `<p:set><p:cBhvr><p:cTn id="${n++}" dur="1" fill="hold"><p:stCondLst><p:cond delay="0"/></p:stCondLst></p:cTn>${tgt(id)}` +
    `<p:attrNameLst><p:attrName>style.visibility</p:attrName></p:attrNameLst></p:cBhvr><p:to><p:strVal val="visible"/></p:to></p:set>`;
  const prop = (id, dur, attr, from, to) =>
    `<p:anim calcmode="lin" valueType="num"><p:cBhvr additive="base"><p:cTn id="${n++}" dur="${dur}" fill="hold"/>${tgt(id)}` +
    `<p:attrNameLst><p:attrName>${attr}</p:attrName></p:attrNameLst></p:cBhvr><p:tavLst>` +
    `<p:tav tm="0"><p:val><p:strVal val="${from}"/></p:val></p:tav><p:tav tm="100000"><p:val><p:strVal val="${to}"/></p:val></p:tav></p:tavLst></p:anim>`;
  const filter = (id, dur, f) =>
    `<p:animEffect transition="in" filter="${f}"><p:cBhvr><p:cTn id="${n++}" dur="${dur}"/>${tgt(id)}</p:cBhvr></p:animEffect>`;

  const PRESET = { fade: [10, 0], fly: [2, 4], wipe: [22, 8], zoom: [53, 16] };
  function effect(kind, id, dur, nodeType) {
    const [presetID, sub] = PRESET[kind];
    const outer = n++;
    let body = setVisible(id);
    if (kind === "fade") body += filter(id, dur, "fade");
    if (kind === "wipe") body += filter(id, dur, "wipe(left)");
    if (kind === "fly") body += prop(id, dur, "ppt_x", "#ppt_x", "#ppt_x") + prop(id, dur, "ppt_y", "#ppt_y+0.08", "#ppt_y") + filter(id, dur, "fade");
    if (kind === "zoom") body += prop(id, dur, "ppt_w", "#ppt_w*0.7", "#ppt_w") + prop(id, dur, "ppt_h", "#ppt_h*0.7", "#ppt_h") + filter(id, dur, "fade");
    return `<p:par><p:cTn id="${outer}" presetID="${presetID}" presetClass="entr" presetSubtype="${sub}" fill="hold" grpId="0" nodeType="${nodeType}">` +
      `<p:stCondLst><p:cond delay="0"/></p:stCondLst><p:childTnLst>${body}</p:childTnLst></p:cTn></p:par>`;
  }

  let groups = "";
  const built = [];
  SEQUENCE.forEach(([start, kind, dur, names], gi) => {
    const gid = n++;
    let effects = "";
    names.forEach((name, k) => {
      const id = ids[name];
      if (!id) throw new Error(`no shape named ${name}`);
      effects += effect(kind, id, dur, gi === 0 && k === 0 ? "afterEffect" : "withEffect");
      built.push(id);
    });
    groups += `<p:par><p:cTn id="${gid}" fill="hold"><p:stCondLst><p:cond delay="${start}"/></p:stCondLst><p:childTnLst>${effects}</p:childTnLst></p:cTn></p:par>`;
  });

  const timing =
    `<p:timing><p:tnLst><p:par><p:cTn id="1" dur="indefinite" restart="never" nodeType="tmRoot"><p:childTnLst>` +
    `<p:seq concurrent="1" nextAc="seek"><p:cTn id="2" dur="indefinite" nodeType="mainSeq"><p:childTnLst>` +
    `<p:par><p:cTn id="${n++}" fill="hold"><p:stCondLst><p:cond delay="indefinite"/><p:cond evt="onBegin" delay="0"><p:tn val="2"/></p:cond></p:stCondLst>` +
    `<p:childTnLst>${groups}</p:childTnLst></p:cTn></p:par>` +
    `</p:childTnLst></p:cTn><p:prevCondLst><p:cond evt="onPrev" delay="0"><p:tgtEl><p:sldTgt/></p:tgtEl></p:cond></p:prevCondLst>` +
    `<p:nextCondLst><p:cond evt="onNext" delay="0"><p:tgtEl><p:sldTgt/></p:tgtEl></p:cond></p:nextCondLst></p:seq>` +
    `</p:childTnLst></p:cTn></p:par></p:tnLst><p:bldLst>` +
    [...new Set(built)].filter((id) => new RegExp(`<p:cNvPr id="${id}"[^>]*>(?:(?!<p:cNvPr).)*?</p:nvSpPr>`, "s").test(xml))
      .map((id) => `<p:bldP spid="${id}" grpId="0" animBg="1"/>`).join("") +
    `</p:bldLst></p:timing>`;
  const transition = `<p:transition spd="med"><p:fade/></p:transition>`;

  xml = xml.replace(/<p:transition[\s\S]*?(<\/p:transition>|\/>)/, "").replace(/<p:timing>[\s\S]*?<\/p:timing>/, "");
  xml = xml.replace("</p:clrMapOvr>", "</p:clrMapOvr>" + transition + timing);
  if (!xml.includes("<p:timing>")) throw new Error("could not place <p:timing>");
  zip.file(slidePath, xml);
  fs.writeFileSync(file, await zip.generateAsync({ type: "nodebuffer", compression: "DEFLATE" }));
}

main().catch((e) => { console.error(e); process.exit(1); });
