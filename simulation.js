// simulation.js - Chemical equilibrium engine, graphing, and molecular visualization

// HTML UI Elements
const acidConcSlider = document.getElementById('acidConcSlider');
const acidConcVal = document.getElementById('acidConcVal');
const baseConcSlider = document.getElementById('baseConcSlider');
const baseConcVal = document.getElementById('baseConcVal');
const pkaSlider = document.getElementById('pkaSlider');
const pkaVal = document.getElementById('pkaVal');
const pkbSlider = document.getElementById('pkbSlider');
const pkbVal = document.getElementById('pkbVal');

const pkaGroup = document.getElementById('pkaGroup');
const pkbGroup = document.getElementById('pkbGroup');

const phDisplay = document.getElementById('phDisplay');
const pohDisplay = document.getElementById('pohDisplay');
const hConcDisplay = document.getElementById('hConcDisplay');
const ohConcDisplay = document.getElementById('ohConcDisplay');
const condDisplay = document.getElementById('condDisplay');

const buretLiquid = document.getElementById('buretLiquid');
const beakerLiquid = document.getElementById('beakerLiquid');
const stopcockValve = document.getElementById('stopcockValve');
const dripLine = document.getElementById('dripLine');
const bulbGlow = document.getElementById('bulbGlow');

const analyticsCanvas = document.getElementById('analyticsChart');
const molecularCanvas = document.getElementById('molecularCanvas');

// Species Legend Labels
const legHA = document.getElementById('legHA');
const legHAText = document.getElementById('legHAText');
const legA = document.getElementById('legA');
const legAText = document.getElementById('legAText');
const legB = document.getElementById('legB');
const legBText = document.getElementById('legBText');
const legBH = document.getElementById('legBH');
const legBHText = document.getElementById('legBHText');

// Simulation Constants
const Kw = 1e-14;
const ACID_VOL = 25.0; // Fixed analyte acid volume in mL
const MAX_BURET_VOL = 50.0; // Buret capacity in mL

// Molecular species radii and colors
const PARTICLE_STYLES = {
  'H': { r: 5, color: '#ff2a70', label: 'H₃O⁺', glow: 'rgba(255, 42, 112, 0.7)' },
  'OH': { r: 5, color: '#3b82f6', label: 'OH⁻', glow: 'rgba(59, 130, 246, 0.7)' },
  'HA': { r: 8, color: '#64748b', label: 'HA', glow: '' },
  'A': { r: 8, color: '#e2e8f0', label: 'A⁻', glow: 'rgba(226, 232, 240, 0.4)' },
  'B': { r: 8, color: '#f59e0b', label: 'B', glow: '' },
  'BH': { r: 8, color: '#a78bfa', label: 'BH⁺', glow: 'rgba(167, 139, 250, 0.4)' }
};

// State Variables
let currentSystem = 'SA-SB'; // SA-SB, SA-WB, WA-SB, WA-WB
let selectedIndicator = 'None'; // None, Phenolphthalein, MethylOrange, Bromothymol
let chartType = 'pH'; // pH, Conductivity
let addedBaseVol = 0.0; // Current base added in mL
let isAutoTitrating = false;
let autoTitrateInterval = null;

// Precalculated curve data for background graphing
let precalcPHCurve = [];
let precalcCondCurve = [];
let activePoints = []; // Titrated points so far

// Particle state lists
let particles = [];
let sparks = [];

// Unified chemical parameters
let chemParams = {
  Ca0: 0.10, // Initial acid concentration
  Cb0: 0.10, // Initial base concentration
  Ka: 1e7,   // Acid dissociation constant (10^7 for SA)
  Kb: 1e7    // Base dissociation constant (10^7 for SB)
};

// -------------------------------------------------------------
// 1. Unified Chemical Equilibrium Solver
// -------------------------------------------------------------

/**
 * Solve for [H+] (x) using charge balance and bisection.
 * x + [BH+] - [OH-] - [A-] = 0
 * where:
 *   [OH-] = Kw / x
 *   [A-] = Ca * Ka / (x + Ka)
 *   [BH+] = Cb * x / (x + Kaconj)    (Kaconj = Kw / Kb)
 */
function calculateEquilibriumState(Ca, Cb, Ka, Kb) {
  const Kaconj = Kw / Kb;
  
  let lowPH = 0;
  let highPH = 14;
  let pH = 7.0;
  
  for (let iter = 0; iter < 50; iter++) {
    pH = (lowPH + highPH) / 2;
    const x = Math.pow(10, -pH); // [H+]
    
    const bh = (Cb * x) / (x + Kaconj);
    const oh = Kw / x;
    const a = (Ca * Ka) / (x + Ka);
    
    const val = x + bh - oh - a;
    
    if (val > 0) {
      lowPH = pH; // x is too high (pH too low)
    } else {
      highPH = pH; // x is too low (pH too high)
    }
  }
  
  const finalPH = pH;
  const H = Math.pow(10, -finalPH);
  const OH = Kw / H;
  
  // Calculate equilibrium species concentrations
  const A = (Ca * Ka) / (H + Ka);
  const HA = Ca - A;
  const BH = (Cb * H) / (H + Kaconj);
  const B = Cb - BH;
  
  // Calculate conductivity relative to ion concentration and mobilities
  // Mobilities in standard units: H+=349.8, OH-=198.6, Na+=50.1, Cl-=76.3, CH3COO-=40.9, NH4+=73.5
  const mobH = 349.8;
  const mobOH = 198.6;
  
  // Assign anion and cation mobilities based on system type
  const mobA = (Ka > 1e3) ? 76.3 : 40.9;     // Cl- (76.3) vs CH3COO- (40.9)
  const mobBH = (Kb > 1e3) ? 50.1 : 73.5;   // Na+ (50.1) vs NH4+ (73.5)
  
  const rawConductivity = (H * mobH) + (OH * mobOH) + (A * mobA) + (BH * mobBH);
  
  return {
    pH: finalPH,
    pOH: 14 - finalPH,
    H: H,
    OH: OH,
    HA: Math.max(0, HA),
    A: Math.max(0, A),
    B: Math.max(0, B),
    BH: Math.max(0, BH),
    conductivity: rawConductivity
  };
}

/**
 * Calculates current concentrations in the beaker after dilution and reaction
 */
function getDilutedConcentrations(vBase) {
  const Vt = ACID_VOL + vBase;
  const Ca = chemParams.Ca0 * (ACID_VOL / Vt);
  const Cb = chemParams.Cb0 * (vBase / Vt);
  return calculateEquilibriumState(Ca, Cb, chemParams.Ka, chemParams.Kb);
}

// -------------------------------------------------------------
// 2. Precalculate Titration Curves
// -------------------------------------------------------------
function precalculateCurves() {
  precalcPHCurve = [];
  precalcCondCurve = [];
  
  // Precalculate values for added base volume from 0 to 100 mL
  const maxV = getMaxVolume();
  const step = maxV / 200;
  
  let refMaxCond = 1.0;
  
  for (let v = 0; v <= maxV + 0.001; v += step) {
    const state = getDilutedConcentrations(v);
    
    // Find maximum conductivity to use as reference scale
    if (v === 0 && currentSystem.startsWith('SA')) {
      refMaxCond = state.conductivity; // Initial strong acid is highest conductivity reference
    } else if (v === 0 && state.conductivity > refMaxCond) {
      refMaxCond = state.conductivity;
    }
    
    precalcPHCurve.push({ volume: v, pH: state.pH });
    precalcCondCurve.push({ volume: v, cond: state.conductivity });
  }
  
  // If weak acid, reference conductivity is scaled by a factor to make it look prominent
  if (!currentSystem.startsWith('SA')) {
    refMaxCond = 15.0; // Standardize conductivity scale for weak acid
  }
  
  // Convert conductivity to relative percentage
  precalcCondCurve.forEach(pt => {
    pt.condPct = Math.min(100, Math.max(0, (pt.cond / refMaxCond) * 100));
  });
}

function getMaxVolume() {
  // Expected equivalence volume: V_eq = V_acid * Ca / Cb
  const Veq = ACID_VOL * chemParams.Ca0 / chemParams.Cb0;
  return Veq * 2; // Graph shows up to 2x equivalence volume
}

// -------------------------------------------------------------
// 3. UI Controls Updates & Custom Graph Drawing
// -------------------------------------------------------------

function setChartType(type) {
  chartType = type;
  document.getElementById('btnShowPH').classList.toggle('active', type === 'pH');
  document.getElementById('btnShowCond').classList.toggle('active', type === 'Conductivity');
  drawGraph();
}

function selectIndicator(indName) {
  selectedIndicator = indName;
  document.querySelectorAll('.ind-btn').forEach(btn => {
    btn.classList.remove('active');
  });
  
  if (indName === 'None') document.getElementById('indNone').classList.add('active');
  else if (indName === 'Phenolphthalein') document.getElementById('indPhth').classList.add('active');
  else if (indName === 'MethylOrange') document.getElementById('indMO').classList.add('active');
  else if (indName === 'Bromothymol') document.getElementById('indBTB').classList.add('active');
  
  updateChemicalState();
}

function updateChemicalState() {
  const state = getDilutedConcentrations(addedBaseVol);
  
  // 1. Update Digital Displays
  phDisplay.textContent = state.pH.toFixed(2);
  pohDisplay.textContent = state.pOH.toFixed(2);
  
  // Format scientific notation for concentrations
  hConcDisplay.innerHTML = formatScientific(state.H);
  ohConcDisplay.innerHTML = formatScientific(state.OH);
  
  // Color the digital pH display based on acidic/neutral/basic
  let phColor = 'var(--color-neutral)';
  let phGlow = 'rgba(0, 235, 208, 0.2)';
  if (state.pH < 6.0) {
    phColor = 'var(--color-acid)';
    phGlow = 'rgba(255, 42, 112, 0.2)';
  } else if (state.pH > 8.0) {
    phColor = 'var(--color-base)';
    phGlow = 'rgba(139, 92, 246, 0.2)';
  }
  phDisplay.style.color = phColor;
  phDisplay.style.textShadow = `0 0 20px ${phGlow}`;
  
  // 2. Update Beaker Color based on indicator and pH
  const beakerColor = getBeakerLiquidColor(state.pH);
  beakerLiquid.style.fill = beakerColor;
  
  // Update liquid level in beaker (starts at 35% height, increases to 85% at max base)
  const maxV = getMaxVolume();
  const fillRatio = addedBaseVol / maxV;
  // SVG coordinate path d: Y starts at 55 (empty) down to 20 (full)
  const yTop = 55 - (35 * fillRatio);
  beakerLiquid.setAttribute('d', `M21 ${yTop} L21 85 A4 4 0 0 0 25 89 L75 89 A4 4 0 0 0 79 85 L79 ${yTop} Z`);
  
  // 3. Update Buret Liquid height
  const buretHeightRatio = (MAX_BURET_VOL - addedBaseVol) / MAX_BURET_VOL;
  // Buret liquid starts at Y = 20, height = 159. When empty, height -> 0, Y -> 179
  const liquidHeight = Math.max(0, 159 * buretHeightRatio);
  const liquidY = 20 + (159 - liquidHeight);
  buretLiquid.setAttribute('y', liquidY);
  buretLiquid.setAttribute('height', liquidHeight);
  
  // 4. Update Conductivity lightbulb brightness
  let relCond = 0;
  if (precalcCondCurve.length > 0) {
    // Map current addedBaseVol to precalcCondCurve percentage
    const step = maxV / 200;
    const idx = Math.min(precalcCondCurve.length - 1, Math.round(addedBaseVol / step));
    relCond = precalcCondCurve[idx] ? precalcCondCurve[idx].condPct : 0;
  }
  condDisplay.textContent = `${Math.round(relCond)}%`;
  
  if (relCond > 1.5) {
    bulbGlow.classList.add('glowing');
    bulbGlow.style.opacity = (relCond / 100) * 0.95 + 0.05;
    bulbGlow.style.filter = `blur(${Math.max(2, 6 * (relCond / 100))}px)`;
  } else {
    bulbGlow.classList.remove('glowing');
    bulbGlow.style.opacity = 0;
  }
  
  // 5. Update graph and active tracer points
  updateActiveTracerPoints();
  drawGraph();
  
  // 6. Update target counts for molecular simulation
  updateMolecularTargets(state);
}

function formatScientific(val) {
  if (val >= 0.01) return val.toFixed(3);
  const exp = Math.floor(Math.log10(val));
  const base = val / Math.pow(10, exp);
  return `${base.toFixed(2)} × 10<sup>${exp}</sup>`;
}

function getBeakerLiquidColor(pH) {
  if (selectedIndicator === 'None') {
    // Neutral plain water blueish tint
    return 'rgba(0, 235, 208, 0.15)';
  }
  
  if (selectedIndicator === 'Phenolphthalein') {
    // Colorless below 8.2, pink above 10
    if (pH < 8.2) return 'rgba(255, 255, 255, 0.05)';
    if (pH > 10.0) return 'rgba(255, 0, 127, 0.65)';
    const pct = (pH - 8.2) / 1.8;
    return `rgba(255, 0, 127, ${pct * 0.65})`;
  }
  
  if (selectedIndicator === 'MethylOrange') {
    // Red below 3.1, Yellow above 4.4
    if (pH < 3.1) return 'rgba(239, 68, 68, 0.6)'; // Red
    if (pH > 4.4) return 'rgba(245, 158, 11, 0.6)'; // Yellow-orange
    // Interpolate hue from red (0) to yellow (45)
    const pct = (pH - 3.1) / 1.3;
    const hue = Math.round(pct * 45);
    return `hsla(${hue}, 100%, 50%, 0.6)`;
  }
  
  if (selectedIndicator === 'Bromothymol') {
    // Yellow below 6.0, Blue above 7.6, Green in buffer
    if (pH < 6.0) return 'rgba(245, 158, 11, 0.6)'; // Yellow
    if (pH > 7.6) return 'rgba(59, 130, 246, 0.6)'; // Blue
    // Interpolate hue from yellow (45) to blue (210)
    const pct = (pH - 6.0) / 1.6;
    const hue = Math.round(45 + pct * 165);
    return `hsla(${hue}, 80%, 45%, 0.6)`;
  }
  
  return 'rgba(0, 235, 208, 0.15)';
}

function updateActiveTracerPoints() {
  activePoints = [];
  const maxV = getMaxVolume();
  
  if (chartType === 'pH') {
    precalcPHCurve.forEach(pt => {
      if (pt.volume <= addedBaseVol + 0.001) {
        activePoints.push({ x: pt.volume, y: pt.pH });
      }
    });
  } else {
    precalcCondCurve.forEach(pt => {
      if (pt.volume <= addedBaseVol + 0.001) {
        activePoints.push({ x: pt.volume, y: pt.condPct });
      }
    });
  }
}

/**
 * Draw custom charts inside the canvas card
 */
function drawGraph() {
  const ctx = analyticsCanvas.getContext('2d');
  
  // Resize canvas to parent bounds
  analyticsCanvas.width = analyticsCanvas.parentElement.clientWidth;
  analyticsCanvas.height = analyticsCanvas.parentElement.clientHeight;
  
  const w = analyticsCanvas.width;
  const h = analyticsCanvas.height;
  
  ctx.clearRect(0, 0, w, h);
  
  const padding = { top: 20, right: 30, bottom: 45, left: 45 };
  const graphWidth = w - padding.left - padding.right;
  const graphHeight = h - padding.top - padding.bottom;
  
  const maxVolume = getMaxVolume();
  const maxValY = (chartType === 'pH') ? 14.0 : 100.0;
  
  // Coordinate converters
  const getX = (vol) => padding.left + (vol / maxVolume) * graphWidth;
  const getY = (val) => padding.top + (1 - (val / maxValY)) * graphHeight;
  
  // 1. Shaded Indicator Transitions (for pH only)
  if (chartType === 'pH' && selectedIndicator !== 'None') {
    let lowerY, upperY, blockColor;
    if (selectedIndicator === 'Phenolphthalein') {
      lowerY = getY(10.0);
      upperY = getY(8.2);
      blockColor = 'rgba(255, 0, 127, 0.06)';
    } else if (selectedIndicator === 'MethylOrange') {
      lowerY = getY(4.4);
      upperY = getY(3.1);
      blockColor = 'rgba(245, 158, 11, 0.06)';
    } else if (selectedIndicator === 'Bromothymol') {
      lowerY = getY(7.6);
      upperY = getY(6.0);
      blockColor = 'rgba(16, 185, 129, 0.06)';
    }
    
    if (blockColor) {
      ctx.fillStyle = blockColor;
      ctx.fillRect(padding.left, lowerY, graphWidth, upperY - lowerY);
      
      // Draw indicator band boundaries
      ctx.strokeStyle = 'rgba(255,255,255,0.03)';
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(padding.left, lowerY);
      ctx.lineTo(padding.left + graphWidth, lowerY);
      ctx.moveTo(padding.left, upperY);
      ctx.lineTo(padding.left + graphWidth, upperY);
      ctx.stroke();
    }
  }
  
  // 2. Draw Gridlines
  ctx.strokeStyle = 'rgba(255, 255, 255, 0.05)';
  ctx.lineWidth = 1;
  ctx.fillStyle = 'rgba(255, 255, 255, 0.4)';
  ctx.font = '10px Outfit';
  ctx.textAlign = 'right';
  ctx.textBaseline = 'middle';
  
  // Horizontal gridlines and Y-axis labels
  const yTicks = (chartType === 'pH') ? [0, 2, 4, 6, 8, 10, 12, 14] : [0, 20, 40, 60, 80, 100];
  yTicks.forEach(tick => {
    const y = getY(tick);
    
    // Gridline
    ctx.beginPath();
    ctx.moveTo(padding.left, y);
    ctx.lineTo(padding.left + graphWidth, y);
    ctx.stroke();
    
    // Label
    ctx.fillText(tick + (chartType === 'Conductivity' ? '%' : ''), padding.left - 8, y);
  });
  
  // Vertical gridlines and X-axis labels
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  const xTicksCount = 5;
  for (let i = 0; i <= xTicksCount; i++) {
    const vol = (maxVolume / xTicksCount) * i;
    const x = getX(vol);
    
    // Gridline
    ctx.beginPath();
    ctx.moveTo(x, padding.top);
    ctx.lineTo(x, padding.top + graphHeight);
    ctx.stroke();
    
    // Label
    ctx.fillText(vol.toFixed(0) + ' mL', x, padding.top + graphHeight + 8);
  }
  
  // Axis lines
  ctx.strokeStyle = 'rgba(255, 255, 255, 0.15)';
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  // X axis
  ctx.moveTo(padding.left, padding.top + graphHeight);
  ctx.lineTo(padding.left + graphWidth, padding.top + graphHeight);
  // Y axis
  ctx.moveTo(padding.left, padding.top);
  ctx.lineTo(padding.left, padding.top + graphHeight);
  ctx.stroke();
  
  // Graph Labels
  ctx.save();
  ctx.font = '11px Outfit';
  ctx.fillStyle = 'rgba(255, 255, 255, 0.6)';
  ctx.textAlign = 'center';
  // X axis label
  ctx.fillText('Volume of Base Added (cm³)', padding.left + graphWidth / 2, padding.top + graphHeight + 28);
  // Y axis label
  ctx.translate(12, padding.top + graphHeight / 2);
  ctx.rotate(-Math.PI / 2);
  ctx.fillText(chartType === 'pH' ? 'pH Value' : 'Relative Conductivity (%)', 0, 0);
  ctx.restore();
  
  // 3. Draw Theoretical Reference Curve (gray outline)
  ctx.strokeStyle = 'rgba(255, 255, 255, 0.15)';
  ctx.lineWidth = 1.5;
  ctx.setLineDash([4, 4]);
  ctx.beginPath();
  
  if (chartType === 'pH') {
    precalcPHCurve.forEach((pt, i) => {
      const cx = getX(pt.volume);
      const cy = getY(pt.pH);
      if (i === 0) ctx.moveTo(cx, cy);
      else ctx.lineTo(cx, cy);
    });
  } else {
    precalcCondCurve.forEach((pt, i) => {
      const cx = getX(pt.volume);
      const cy = getY(pt.condPct);
      if (i === 0) ctx.moveTo(cx, cy);
      else ctx.lineTo(cx, cy);
    });
  }
  ctx.stroke();
  ctx.setLineDash([]); // reset line dash
  
  // 4. Draw Active Tracer Curve (Glowing Neon Line)
  if (activePoints.length > 0) {
    ctx.strokeStyle = (chartType === 'pH') ? 'var(--color-acid)' : 'var(--color-conductivity)';
    ctx.lineWidth = 3;
    ctx.shadowBlur = 8;
    ctx.shadowColor = ctx.strokeStyle;
    ctx.beginPath();
    activePoints.forEach((pt, i) => {
      const cx = getX(pt.x);
      const cy = getY(pt.y);
      if (i === 0) ctx.moveTo(cx, cy);
      else ctx.lineTo(cx, cy);
    });
    ctx.stroke();
    ctx.shadowBlur = 0; // reset shadow glow
    
    // Draw current active dot
    const head = activePoints[activePoints.length - 1];
    const headX = getX(head.x);
    const headY = getY(head.y);
    
    ctx.fillStyle = ctx.strokeStyle;
    ctx.beginPath();
    ctx.arc(headX, headY, 5, 0, Math.PI * 2);
    ctx.fill();
    ctx.strokeStyle = '#fff';
    ctx.lineWidth = 1.5;
    ctx.stroke();
  }
  
  // 5. Draw Equivalence & Half-Equivalence Annotations
  const Veq = ACID_VOL * chemParams.Ca0 / chemParams.Cb0;
  
  // Only draw if we have titrated up to that point
  if (addedBaseVol >= Veq - 0.1) {
    const eqState = getDilutedConcentrations(Veq);
    const eqX = getX(Veq);
    const eqY = (chartType === 'pH') ? getY(eqState.pH) : getY(precalcCondCurve[Math.round(precalcCondCurve.length / 2)].condPct);
    
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.3)';
    ctx.lineWidth = 1;
    ctx.setLineDash([2, 2]);
    
    // Vertical equivalence indicator
    ctx.beginPath();
    ctx.moveTo(eqX, padding.top);
    ctx.lineTo(eqX, padding.top + graphHeight);
    ctx.stroke();
    
    ctx.setLineDash([]);
    
    // Equivalence point annotation dot
    ctx.fillStyle = 'var(--color-accent)';
    ctx.beginPath();
    ctx.arc(eqX, eqY, 6, 0, Math.PI * 2);
    ctx.fill();
    ctx.strokeStyle = '#ffffff';
    ctx.lineWidth = 1.5;
    ctx.stroke();
    
    // Text label
    ctx.fillStyle = '#ffffff';
    ctx.font = 'bold 9px Outfit';
    ctx.textAlign = 'left';
    ctx.fillText('Equivalence Pt', eqX + 8, eqY - 4);
    ctx.font = 'normal 9px Outfit';
    ctx.fillStyle = 'var(--text-secondary)';
    
    if (chartType === 'pH') {
      ctx.fillText(`Vol: ${Veq.toFixed(1)} mL, pH: ${eqState.pH.toFixed(2)}`, eqX + 8, eqY + 6);
    }
  }
  
  // For Weak Acid titrations: Half-Equivalence Point (V_eq / 2)
  if (currentSystem.startsWith('WA') && addedBaseVol >= (Veq / 2) - 0.1 && chartType === 'pH') {
    const halfV = Veq / 2;
    const halfState = getDilutedConcentrations(halfV);
    const hX = getX(halfV);
    const hY = getY(halfState.pH);
    
    ctx.fillStyle = 'var(--color-success)';
    ctx.beginPath();
    ctx.arc(hX, hY, 5, 0, Math.PI * 2);
    ctx.fill();
    ctx.strokeStyle = '#ffffff';
    ctx.lineWidth = 1.5;
    ctx.stroke();
    
    ctx.fillStyle = '#ffffff';
    ctx.font = 'bold 9px Outfit';
    ctx.textAlign = 'right';
    ctx.fillText('Half-Equivalence', hX - 8, hY - 4);
    ctx.font = 'italic bold 9px Outfit';
    ctx.fillStyle = 'var(--color-neutral)';
    ctx.fillText(`pH = pKₐ = ${halfState.pH.toFixed(2)}`, hX - 8, hY + 6);
  }
}

// -------------------------------------------------------------
// 4. Molecular-Scale Canvas Simulation
// -------------------------------------------------------------

class Particle {
  constructor(type, x, y) {
    this.type = type;
    this.x = x;
    this.y = y;
    
    const style = PARTICLE_STYLES[type];
    this.r = style.r;
    this.color = style.color;
    this.glow = style.glow;
    
    // Random thermal velocities
    const speed = (type === 'H' || type === 'OH') ? 2.2 : 1.2; // Protons jump faster
    const angle = Math.random() * Math.PI * 2;
    this.vx = Math.cos(angle) * speed;
    this.vy = Math.sin(angle) * speed;
    
    this.id = Math.random();
    this.fade = 1.0;
    this.isFading = false;
  }
  
  draw(ctx) {
    ctx.save();
    ctx.globalAlpha = this.fade;
    
    // Draw halo glow if active
    if (this.glow && this.fade > 0.1) {
      ctx.shadowBlur = 8;
      ctx.shadowColor = this.color;
    }
    
    ctx.beginPath();
    ctx.arc(this.x, this.y, this.r, 0, Math.PI * 2);
    ctx.fillStyle = this.color;
    ctx.fill();
    
    // Draw border
    ctx.strokeStyle = 'rgba(255,255,255,0.2)';
    ctx.lineWidth = 1.0;
    ctx.stroke();
    
    ctx.shadowBlur = 0; // reset
    
    // Draw chemical charge / label text inside particle
    ctx.fillStyle = '#ffffff';
    ctx.font = 'bold 8px Outfit';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    
    let symbol = '';
    if (this.type === 'H') symbol = '+';
    else if (this.type === 'OH') symbol = '-';
    else if (this.type === 'HA') symbol = 'HA';
    else if (this.type === 'A') symbol = 'A⁻';
    else if (this.type === 'B') symbol = 'B';
    else if (this.type === 'BH') symbol = 'BH⁺';
    
    ctx.fillText(symbol, this.x, this.y);
    ctx.restore();
  }
  
  update(w, h) {
    if (this.isFading) {
      this.fade -= 0.05;
      return;
    }
    
    this.x += this.vx;
    this.y += this.vy;
    
    // Boundary bounces
    if (this.x - this.r < 0) {
      this.x = this.r;
      this.vx = -this.vx;
    } else if (this.x + this.r > w) {
      this.x = w - this.r;
      this.vx = -this.vx;
    }
    
    if (this.y - this.r < 0) {
      this.y = this.r;
      this.vy = -this.vy;
    } else if (this.y + this.r > h) {
      this.y = h - this.r;
      this.vy = -this.vy;
    }
  }
}

// Spark effect particles
function createSparks(x, y, color) {
  for (let i = 0; i < 8; i++) {
    const angle = Math.random() * Math.PI * 2;
    const speed = 1.0 + Math.random() * 2.0;
    sparks.push({
      x: x,
      y: y,
      vx: Math.cos(angle) * speed,
      vy: Math.sin(angle) * speed,
      life: 15,
      maxLife: 15,
      color: color || '#00ebd0'
    });
  }
}

function updateAndDrawSparks(ctx) {
  for (let i = sparks.length - 1; i >= 0; i--) {
    const s = sparks[i];
    s.x += s.vx;
    s.y += s.vy;
    s.life--;
    
    if (s.life <= 0) {
      sparks.splice(i, 1);
      continue;
    }
    
    ctx.save();
    ctx.globalAlpha = s.life / s.maxLife;
    ctx.beginPath();
    ctx.arc(s.x, s.y, 2, 0, Math.PI * 2);
    ctx.fillStyle = s.color;
    ctx.shadowBlur = 6;
    ctx.shadowColor = s.color;
    ctx.fill();
    ctx.restore();
  }
}

// Particle target limits computed from chemical equilibrium concentrations
let targetCounts = { H: 0, OH: 0, HA: 0, A: 0, B: 0, BH: 0 };

function updateMolecularTargets(state) {
  const totalConc = state.H + state.OH + state.HA + state.A + state.B + state.BH;
  const N_MAX = 42; // Maximum particles to avoid clutter
  
  if (totalConc > 0) {
    targetCounts.H = Math.round(N_MAX * state.H / totalConc);
    targetCounts.OH = Math.round(N_MAX * state.OH / totalConc);
    targetCounts.HA = Math.round(N_MAX * state.HA / totalConc);
    targetCounts.A = Math.round(N_MAX * state.A / totalConc);
    targetCounts.B = Math.round(N_MAX * state.B / totalConc);
    targetCounts.BH = Math.round(N_MAX * state.BH / totalConc);
    
    // Ensure active types exist visual representation even if low
    if (state.H > 1e-5 && targetCounts.H === 0) targetCounts.H = 1;
    if (state.OH > 1e-5 && targetCounts.OH === 0) targetCounts.OH = 1;
  }
}

function adjustParticlePopulation() {
  const w = molecularCanvas.width;
  
  // Count current active non-fading particles
  const currentCounts = { H: 0, OH: 0, HA: 0, A: 0, B: 0, BH: 0 };
  particles.forEach(p => {
    if (!p.isFading) currentCounts[p.type]++;
  });
  
  // For each type, spawn or remove
  Object.keys(targetCounts).forEach(type => {
    const diff = targetCounts[type] - currentCounts[type];
    
    if (diff > 0) {
      // Spawn new particles from the top center (entering from buret)
      for (let i = 0; i < diff; i++) {
        const x = w / 2 + (Math.random() - 0.5) * 20;
        const y = 10;
        particles.push(new Particle(type, x, y));
      }
    } else if (diff < 0) {
      // Flag oldest particles of this type to fade out
      let removed = 0;
      for (let i = 0; i < particles.length && removed < Math.abs(diff); i++) {
        const p = particles[i];
        if (p.type === type && !p.isFading) {
          p.isFading = true;
          removed++;
        }
      }
    }
  });
}

function checkMolecularCollisions() {
  // Check for reactive neutralization: H3O+ + OH- -> H2O (Both disappear)
  for (let i = 0; i < particles.length; i++) {
    for (let j = i + 1; j < particles.length; j++) {
      const p1 = particles[i];
      const p2 = particles[j];
      
      if (p1.isFading || p2.isFading) continue;
      
      const dx = p2.x - p1.x;
      const dy = p2.y - p1.y;
      const dist = Math.sqrt(dx * dx + dy * dy);
      const minDist = p1.r + p2.r;
      
      if (dist < minDist) {
        // Resolve overlap bouncing
        const overlap = minDist - dist;
        const nx = dx / dist;
        const ny = dy / dist;
        
        p1.x -= nx * overlap * 0.5;
        p1.y -= ny * overlap * 0.5;
        p2.x += nx * overlap * 0.5;
        p2.y += ny * overlap * 0.5;
        
        // Elastic bounce logic
        const kx = p1.vx - p2.vx;
        const ky = p1.vy - p2.vy;
        const vn = kx * nx + ky * ny;
        if (vn > 0) {
          const impulse = vn; // masses assumed equivalent
          p1.vx -= impulse * nx;
          p1.vy -= impulse * ny;
          p2.vx += impulse * nx;
          p2.vy += impulse * ny;
          
          // --- Chemical Reactions upon collision ---
          // 1. Neutralization: H+ (H) + OH- (OH) -> Water
          if ((p1.type === 'H' && p2.type === 'OH') || (p1.type === 'OH' && p2.type === 'H')) {
            const rx = (p1.x + p2.x) / 2;
            const ry = (p1.y + p2.y) / 2;
            createSparks(rx, ry, 'var(--color-neutral)');
            
            // Remove both immediately
            p1.isFading = true;
            p2.isFading = true;
            return;
          }
          
          // 2. Recombination: H+ + A- -> HA (If we need more HA and have excess H & A)
          const isH_A = (p1.type === 'H' && p2.type === 'A') || (p1.type === 'A' && p2.type === 'H');
          if (isH_A && currentCounts('HA') < targetCounts.HA) {
            const rx = (p1.x + p2.x) / 2;
            const ry = (p1.y + p2.y) / 2;
            createSparks(rx, ry, 'var(--color-acid)');
            
            p1.isFading = true;
            p2.isFading = true;
            // Spawn intermediate neutral acid
            particles.push(new Particle('HA', rx, ry));
            return;
          }
          
          // 3. Recombination: OH- + BH+ -> B (If we need more B)
          const isOH_BH = (p1.type === 'OH' && p2.type === 'BH') || (p1.type === 'BH' && p2.type === 'OH');
          if (isOH_BH && currentCounts('B') < targetCounts.B) {
            const rx = (p1.x + p2.x) / 2;
            const ry = (p1.y + p2.y) / 2;
            createSparks(rx, ry, 'var(--color-base)');
            
            p1.isFading = true;
            p2.isFading = true;
            particles.push(new Particle('B', rx, ry));
            return;
          }
        }
      }
    }
  }
}

function currentCounts(type) {
  let c = 0;
  particles.forEach(p => {
    if (p.type === type && !p.isFading) c++;
  });
  return c;
}

/**
 * Main rendering loop for Molecular Canvas
 */
function animateMolecularCanvas() {
  const ctx = molecularCanvas.getContext('2d');
  
  // Match sizes to css box
  molecularCanvas.width = molecularCanvas.parentElement.clientWidth;
  molecularCanvas.height = molecularCanvas.parentElement.clientHeight;
  
  const w = molecularCanvas.width;
  const h = molecularCanvas.height;
  
  ctx.clearRect(0, 0, w, h);
  
  // Adjust particle populations to match chemical concentrations
  adjustParticlePopulation();
  
  // Physics updates and checks
  checkMolecularCollisions();
  
  // Render and update particles
  for (let i = particles.length - 1; i >= 0; i--) {
    const p = particles[i];
    p.update(w, h);
    
    if (p.fade <= 0.0) {
      particles.splice(i, 1);
      continue;
    }
    
    p.draw(ctx);
  }
  
  // Sparks
  updateAndDrawSparks(ctx);
  
  requestAnimationFrame(animateMolecularCanvas);
}

// -------------------------------------------------------------
// 5. Simulation Action Handles
// -------------------------------------------------------------

function loadPreset(system) {
  currentSystem = system;
  addedBaseVol = 0.0;
  if (isAutoTitrating) toggleAutoTitration();
  
  // Reset presets active buttons
  document.querySelectorAll('.preset-btn').forEach(btn => {
    btn.classList.remove('active');
  });
  
  // Update configuration sliders based on selected system type
  if (system === 'SA-SB') {
    document.getElementById('presetSAsb').classList.add('active');
    chemParams.Ka = 1e7;
    chemParams.Kb = 1e7;
    pkaSlider.disabled = true;
    pkbSlider.disabled = true;
    pkaGroup.classList.add('disabled');
    pkbGroup.classList.add('disabled');
    
    // Set typical indicator
    selectIndicator('Phenolphthalein');
    
    // Set Legend labels
    legHA.style.display = 'flex';
    legHAText.textContent = 'HCl (Strong Acid)';
    legA.style.display = 'flex';
    legAText.textContent = 'Cl⁻ (Spectator Anion)';
    legB.style.display = 'none';
    legBH.style.display = 'flex';
    legBHText.textContent = 'Na⁺ (Spectator Cation)';
  } 
  else if (system === 'SA-WB') {
    document.getElementById('presetSAwb').classList.add('active');
    chemParams.Ka = 1e7;
    chemParams.Kb = Math.pow(10, -parseFloat(pkbSlider.value));
    pkaSlider.disabled = true;
    pkbSlider.disabled = false;
    pkaGroup.classList.add('disabled');
    pkbGroup.classList.remove('disabled');
    
    selectIndicator('MethylOrange'); // Acidic equivalence point
    
    legHA.style.display = 'flex';
    legHAText.textContent = 'HCl (Strong Acid)';
    legA.style.display = 'flex';
    legAText.textContent = 'Cl⁻ (Spectator Anion)';
    legB.style.display = 'flex';
    legBText.textContent = 'NH₃ (Weak Base Molecule)';
    legBH.style.display = 'flex';
    legBHText.textContent = 'NH₄⁺ (Conjugate Acid)';
  } 
  else if (system === 'WA-SB') {
    document.getElementById('presetWAsb').classList.add('active');
    chemParams.Ka = Math.pow(10, -parseFloat(pkaSlider.value));
    chemParams.Kb = 1e7;
    pkaSlider.disabled = false;
    pkbSlider.disabled = true;
    pkaGroup.classList.remove('disabled');
    pkbGroup.classList.add('disabled');
    
    selectIndicator('Phenolphthalein'); // Basic equivalence point
    
    legHA.style.display = 'flex';
    legHAText.textContent = 'CH₃COOH (Weak Acid Molecule)';
    legA.style.display = 'flex';
    legAText.textContent = 'CH₃COO⁻ (Conjugate Base)';
    legB.style.display = 'none';
    legBH.style.display = 'flex';
    legBHText.textContent = 'Na⁺ (Spectator Cation)';
  } 
  else if (system === 'WA-WB') {
    document.getElementById('presetWAwb').classList.add('active');
    chemParams.Ka = Math.pow(10, -parseFloat(pkaSlider.value));
    chemParams.Kb = Math.pow(10, -parseFloat(pkbSlider.value));
    pkaSlider.disabled = false;
    pkbSlider.disabled = false;
    pkaGroup.classList.remove('disabled');
    pkbGroup.classList.remove('disabled');
    
    selectIndicator('None'); // Weak-weak has no sharp equivalence jump, indicator useless
    
    legHA.style.display = 'flex';
    legHAText.textContent = 'CH₃COOH (Weak Acid Molecule)';
    legA.style.display = 'flex';
    legAText.textContent = 'CH₃COO⁻ (Conjugate Base)';
    legB.style.display = 'flex';
    legBText.textContent = 'NH₃ (Weak Base Molecule)';
    legBH.style.display = 'flex';
    legBHText.textContent = 'NH₄⁺ (Conjugate Acid)';
  }
  
  particles = []; // Clear visual particles
  precalculateCurves();
  updateChemicalState();
}

function updateSliders() {
  chemParams.Ca0 = parseFloat(acidConcSlider.value);
  chemParams.Cb0 = parseFloat(baseConcSlider.value);
  acidConcVal.textContent = chemParams.Ca0.toFixed(2) + ' M';
  baseConcVal.textContent = chemParams.Cb0.toFixed(2) + ' M';
  
  if (currentSystem.startsWith('WA')) {
    chemParams.Ka = Math.pow(10, -parseFloat(pkaSlider.value));
    pkaVal.textContent = parseFloat(pkaSlider.value).toFixed(2);
  }
  if (currentSystem.endsWith('WB')) {
    chemParams.Kb = Math.pow(10, -parseFloat(pkbSlider.value));
    pkbVal.textContent = parseFloat(pkbSlider.value).toFixed(2);
  }
  
  addedBaseVol = 0.0;
  if (isAutoTitrating) toggleAutoTitration();
  precalculateCurves();
  updateChemicalState();
}

function addDrop() {
  const maxV = getMaxVolume();
  if (addedBaseVol >= maxV) return;
  
  addedBaseVol = Math.min(maxV, addedBaseVol + 0.5);
  triggerDripAnimation();
  updateChemicalState();
}

function addBulk() {
  const maxV = getMaxVolume();
  if (addedBaseVol >= maxV) return;
  
  addedBaseVol = Math.min(maxV, addedBaseVol + 5.0);
  triggerDripAnimation();
  updateChemicalState();
}

function triggerDripAnimation() {
  // Rotate stopcock valve indicator to open (90deg) and animate drip
  stopcockValve.style.transform = 'rotate(90deg)';
  stopcockValve.style.transformOrigin = '20px 188px';
  dripLine.classList.add('flowing');
  
  setTimeout(() => {
    stopcockValve.style.transform = 'rotate(0deg)';
    dripLine.classList.remove('flowing');
  }, 400);
}

function toggleAutoTitration() {
  const btn = document.getElementById('autoTitrateBtn');
  const maxV = getMaxVolume();
  
  if (isAutoTitrating) {
    // Pause
    isAutoTitrating = false;
    clearInterval(autoTitrateInterval);
    btn.innerHTML = '<i class="fa-solid fa-play"></i> Auto-Titrate';
    stopcockValve.style.transform = 'rotate(0deg)';
    dripLine.classList.remove('flowing');
  } else {
    // Start
    if (addedBaseVol >= maxV) addedBaseVol = 0.0; // auto-reset at bounds
    isAutoTitrating = true;
    btn.innerHTML = '<i class="fa-solid fa-pause"></i> Pause';
    stopcockValve.style.transform = 'rotate(90deg)';
    stopcockValve.style.transformOrigin = '20px 188px';
    dripLine.classList.add('flowing');
    
    autoTitrateInterval = setInterval(() => {
      if (addedBaseVol >= maxV) {
        toggleAutoTitration(); // stop
        return;
      }
      addedBaseVol = Math.min(maxV, addedBaseVol + 0.5);
      updateChemicalState();
    }, 150);
  }
}

function resetTitration() {
  addedBaseVol = 0.0;
  if (isAutoTitrating) toggleAutoTitration();
  particles = [];
  updateChemicalState();
}

// -------------------------------------------------------------
// 6. Hook Event Listeners and Initialize
// -------------------------------------------------------------

window.addEventListener('load', () => {
  // Attach event listeners
  acidConcSlider.addEventListener('input', updateSliders);
  baseConcSlider.addEventListener('input', updateSliders);
  pkaSlider.addEventListener('input', updateSliders);
  pkbSlider.addEventListener('input', updateSliders);
  
  // Initialize preset
  loadPreset('SA-SB');
  
  // Start animation loop for molecular view
  animateMolecularCanvas();
});

window.addEventListener('resize', () => {
  drawGraph();
});
