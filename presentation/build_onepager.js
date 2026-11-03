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
const { FaGasPump, FaShieldAlt, FaExclamationTriangle, FaBalanceScale, FaDraftingCompass, FaCogs, FaCubes } = require("react-icons/fa");

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
      { placeholder: { options: { name: "title", type: "title", x: 0.5, y: 0.72, w: 6.1, h: 1.2,
        fontFace: THEME.headFontFace, fontSize: 34, bold: true, color: C.background1, valign: "top", align: "left", margin: 0 },
        text: "Title" } },
      { text: { text: "github.com/SimoneConiglio/Atkinson-cycle-engine-optimisation-",
        options: { x: 6.83, y: 7.02, w: 6.0, h: 0.25, fontSize: 10, color: C.background2, align: "right", margin: 0 } } },
    ],
  });

  pres.addSection({ title: "Overview" });
  const s = pres.addSlide({ masterName: "ONE_PAGER_DARK", sectionTitle: "Overview" });
  s.transition = undefined;

  // ---- header ----
  s.addText("SHELL ECO-MARATHON  ·  MULTIDISCIPLINARY DESIGN", {
    objectName: "kicker", isTextBox: true, x: 0.5, y: 0.4, w: 6.1, h: 0.3, margin: 0,
    fontSize: 12, bold: true, charSpacing: 2, color: C.accent1 });
  s.addText("Range-optimal Atkinson engine linkage", { placeholder: "title" });

  const label = (name, text, x, y, w) => s.addText(text, { objectName: name, isTextBox: true, x, y, w, h: 0.28, margin: 0,
    fontSize: 12, bold: true, charSpacing: 2, color: C.accent1 });

  // ---- the problem ----
  label("problem-label", "THE PROBLEM", 0.5, 2.0, 6.1);
  s.addText([
    { text: "Size an extended-expansion linkage — eleven dimensions and a gear pair, one piston with two top dead centres and two different strokes — so a car goes " },
    { text: "as far as possible on a fixed quantity of fuel", options: { bold: true, color: C.background1 } },
    { text: ". Thermodynamics wants a long expansion stroke, friction wants short lever arms, inertia couples masses to loads, and the car pays for every kilogram." },
  ], { objectName: "problem-text", isTextBox: true, x: 0.5, y: 2.3, w: 6.1, h: 1.22, margin: 0,
    fontSize: 14, color: C.background2, valign: "top" });

  // ---- highlights ----
  label("findings-label", "HIGHLIGHTS", 0.5, 3.62, 6.1);
  const iconColor = THEME.colors.dk1;
  const findings = [
    { ic: FaExclamationTriangle, c: C.accent3, t: "The quasi-static optimum is unbuildable above 2000 rpm" },
    { ic: FaBalanceScale, c: C.accent1, t: "The topology is worth +17.6 %, mostly friction and flywheel" },
    { ic: FaShieldAlt, c: C.accent2, t: "Reliability-based design: 3395 km/L at Pf = 10⁻³" },
  ];
  for (let i = 0; i < findings.length; i++) {
    const f = findings[i], y = 3.94 + i * 0.38;
    s.addShape(pres.shapes.OVAL, { objectName: `finding-${i}-circle`, x: 0.5, y, w: 0.3, h: 0.3,
      fill: { color: f.c }, line: { color: f.c, width: 0 } });
    s.addImage({ objectName: `finding-${i}-icon`, data: await icon(f.ic, iconColor), x: 0.575, y: y + 0.075, w: 0.15, h: 0.15, altText: f.t });
    s.addText(f.t, { objectName: `finding-${i}-text`, isTextBox: true, x: 0.95, y, w: 5.65, h: 0.3, margin: 0,
      fontSize: 14, color: C.background1, valign: "middle" });
  }

  // ---- animated mechanism ----
  const gw = 5.6, gh = 3.5, gx = 12.83 - gw, gy = 0.4;
  s.addShape(pres.shapes.ROUNDED_RECTANGLE, { objectName: "gif-frame", x: gx, y: gy, w: gw, h: gh,
    rectRadius: 0.12, fill: { color: C.background1 }, line: { color: C.background1, width: 0 } });
  s.addImage({ objectName: "mechanism-gif", path: path.join(ROOT, "docs/figures/exlink.gif"),
    x: gx + 0.06, y: gy + 0.0375, w: gw - 0.12, h: (gw - 0.12) * 0.625,
    altText: "The optimised EX-link linkage turning through one cycle, with piston motion, pressure-volume cycle and crankshaft torque" });

  // ---- stat tiles (one row under the GIF) ----
  const stats = [
    { v: "3395", u: "km/L", l: "range", c: C.accent1 },
    { v: "+17.6", u: "%", l: "vs conventional", c: C.accent2 },
    { v: "+19.3", u: "%", l: "over 4 speeds", c: C.accent2 },
    { v: "12.9", u: "kg", l: "engine (vs 16.9)", c: C.accent1 },
  ];
  const sx0 = 6.85, sw = (12.83 - sx0 - 3 * 0.15) / 4, sy = 4.05, sh = 0.8;
  stats.forEach((st, i) => {
    const x = sx0 + i * (sw + 0.15);
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { objectName: `stat-${i}-bg`, x, y: sy, w: sw, h: sh,
      rectRadius: 0.08, fill: { color: C.text2 }, line: { color: C.text2, width: 0 } });
    s.addText([
      { text: st.v, options: { fontSize: 22, bold: true, color: st.c, fontFace: THEME.headFontFace } },
      { text: " " + st.u, options: { fontSize: 11, bold: true, color: st.c } },
    ], { objectName: `stat-${i}-value`, isTextBox: true, x: x + 0.12, y: sy + 0.06, w: sw - 0.2, h: 0.4, margin: 0, valign: "middle" });
    s.addText(st.l, { objectName: `stat-${i}-label`, isTextBox: true, x: x + 0.12, y: sy + 0.47, w: sw - 0.2, h: 0.24,
      margin: 0, fontSize: 10, color: C.background2, valign: "top" });
  });

  // ---- the disciplines ----
  label("disc-label", "THE DISCIPLINES  ·  FIVE COUPLED MODELS IN ONE GEMSEO MDF SCENARIO", 0.5, 5.2, 9.0);
  s.addText("teal outline = solved together as a fixed point (Gauss–Seidel MDA)", {
    objectName: "disc-legend", isTextBox: true, x: 8.33, y: 5.2, w: 4.5, h: 0.28, margin: 0,
    fontSize: 10, italic: true, color: C.background2, align: "right" });
  const disc = [
    { ic: FaDraftingCompass, c: C.accent1, h: "Kinematics & cycle", b: "strokes, two TDCs, Atkinson cycle, indicated efficiency" },
    { ic: FaCogs, c: C.accent2, h: "Dynamics", b: "gas + inertia loads at speed: members, torque, bearings", loop: true },
    { ic: FaCubes, c: C.accent2, h: "Structure", b: "yield, fatigue, buckling set the sections and the mass", loop: true },
    { ic: FaGasPump, c: C.accent1, h: "Range", b: "friction, mass budget, gear pair, burn-and-coast car → km/L" },
    { ic: FaShieldAlt, c: C.accent3, h: "Reliability", b: "FORM + sampling over 17 uncertain parameters → Pf" },
  ];
  const gap = 0.36, dw = (12.33 - 4 * gap) / 5, dy = 5.52, dh = 1.35;
  for (let i = 0; i < disc.length; i++) {
    const d = disc[i], x = 0.5 + i * (dw + gap);
    s.addShape(pres.shapes.ROUNDED_RECTANGLE, { objectName: `disc-${i}-bg`, x, y: dy, w: dw, h: dh, rectRadius: 0.08,
      fill: { color: C.text2 }, line: d.loop ? { color: C.accent2, width: 1.5 } : { color: C.text2, width: 0 } });
    s.addShape(pres.shapes.OVAL, { objectName: `disc-${i}-circle`, x: x + 0.15, y: dy + 0.15, w: 0.38, h: 0.38,
      fill: { color: d.c }, line: { color: d.c, width: 0 } });
    s.addImage({ objectName: `disc-${i}-icon`, data: await icon(d.ic, iconColor), x: x + 0.245, y: dy + 0.245, w: 0.19, h: 0.19, altText: d.h });
    s.addText(d.h, { objectName: `disc-${i}-head`, isTextBox: true, x: x + 0.63, y: dy + 0.15, w: dw - 0.73, h: 0.38, margin: 0,
      fontSize: 13, bold: true, color: C.background1, valign: "middle" });
    s.addText(d.b, { objectName: `disc-${i}-body`, isTextBox: true, x: x + 0.15, y: dy + 0.64, w: dw - 0.27, h: 0.66, margin: 0,
      fontSize: 11, color: C.background2, valign: "top" });
    if (i < disc.length - 1) {
      s.addText(i === 1 ? "⇄" : "→", { objectName: `disc-${i}-arrow`, isTextBox: true, x: x + dw, y: dy + 0.45, w: gap, h: 0.45,
        margin: 0, fontSize: 20, bold: true, align: "center", valign: "middle", color: i === 1 ? C.accent2 : C.background2 });
    }
  }

  s.addNotes("The problem: size an extended-expansion (EX-link) linkage, with eleven dimensions and a gear pair, so a Shell Eco-marathon car goes as far as possible on a fixed quantity of fuel. Four effects conflict: thermodynamics wants a long expansion stroke, friction wants short lever arms, inertia couples part masses to loads, and the car carries every kilogram. Five disciplines are coupled in one GEMSEO MDF scenario: kinematics and cycle, dynamics, structure (dynamics and structure solved together by a Gauss-Seidel MDA because loads set sections and sections set inertia), range, and reliability. Highlights: the quasi-static optimum is unbuildable above 2000 rpm; the topology is worth 17.6 percent over a conventional engine, mostly through friction and flywheel mass; and reliability-based design gives 3395 km/L at a failure probability of 1e-3.");

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
  [900,  "zoom", 800, ["gif-frame", "mechanism-gif"]],
  [1200, "fade", 500, ["problem-label"]],
  [1400, "fade", 900, ["problem-text"]],
  [2600, "fade", 500, ["disc-label", "disc-legend"]],
  [2900, "fly",  500, ["disc-0-bg", "disc-0-circle", "disc-0-icon", "disc-0-head", "disc-0-body", "disc-0-arrow"]],
  [3250, "fly",  500, ["disc-1-bg", "disc-1-circle", "disc-1-icon", "disc-1-head", "disc-1-body", "disc-1-arrow"]],
  [3600, "fly",  500, ["disc-2-bg", "disc-2-circle", "disc-2-icon", "disc-2-head", "disc-2-body", "disc-2-arrow"]],
  [3950, "fly",  500, ["disc-3-bg", "disc-3-circle", "disc-3-icon", "disc-3-head", "disc-3-body", "disc-3-arrow"]],
  [4300, "fly",  500, ["disc-4-bg", "disc-4-circle", "disc-4-icon", "disc-4-head", "disc-4-body"]],
  [4850, "fade", 500, ["findings-label"]],
  [5100, "wipe", 500, ["finding-0-circle", "finding-0-icon", "finding-0-text"]],
  [5400, "wipe", 500, ["finding-1-circle", "finding-1-icon", "finding-1-text"]],
  [5700, "wipe", 500, ["finding-2-circle", "finding-2-icon", "finding-2-text"]],
  [6100, "fly",  450, ["stat-0-bg", "stat-0-value", "stat-0-label"]],
  [6300, "fly",  450, ["stat-1-bg", "stat-1-value", "stat-1-label"]],
  [6500, "fly",  450, ["stat-2-bg", "stat-2-value", "stat-2-label"]],
  [6700, "fly",  450, ["stat-3-bg", "stat-3-value", "stat-3-label"]],
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
