/* =====================================================================
   STATISTICS UTILITIES
   ===================================================================== */

// Error function — Abramowitz & Stegun approximation (max error 1.5×10⁻⁷)
function erf(x) {
  const sign = x < 0 ? -1 : 1;
  x = Math.abs(x);
  const t = 1 / (1 + 0.3275911 * x);
  const p = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
  return sign * (1 - p * Math.exp(-x * x));
}

// Standard normal CDF: Φ(z) = P(Z ≤ z)
function phi(z) {
  return 0.5 * (1 + erf(z / Math.SQRT2));
}

// Inverse standard normal CDF — Peter Acklam's rational approximation
function phiInv(p) {
  const a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
              1.383577518672690e+02, -3.066479806614716e+01,  2.506628277459239e+00];
  const b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
              6.680131188771972e+01, -1.328068155288572e+01];
  const c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
             -2.549732539343734e+00,  4.374664141464968e+00,  2.938163982698783e+00];
  const d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00];

  if (p <= 0) return -Infinity;
  if (p >= 1) return  Infinity;

  const pLow = 0.02425;
  const pHigh = 1 - pLow;
  let q, r;

  if (p < pLow) {
    q = Math.sqrt(-2 * Math.log(p));
    return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
           ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  } else if (p <= pHigh) {
    q = p - 0.5; r = q * q;
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
  } else {
    q = Math.sqrt(-2 * Math.log(1 - p));
    return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
             ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  }
}

// Pre-computed one-sided critical values
const ZCRIT = {
  95:  phiInv(0.95),  // ≈ 1.6449
  98:  phiInv(0.98),  // ≈ 2.0537
  99:  phiInv(0.99),  // ≈ 2.3263
  p80: phiInv(0.80),  // ≈ 0.8416
  p90: phiInv(0.90),  // ≈ 1.2816
};

// Normal probability density function
function normalPDF(z) {
  return Math.exp(-0.5 * z * z) / Math.sqrt(2 * Math.PI);
}

/* =====================================================================
   HELPERS
   ===================================================================== */

function el(id)         { return document.getElementById(id); }
function setText(id, v) { el(id).textContent = v; }
function setHTML(id, v) { el(id).innerHTML   = v; }

function fmtPct(v, dp = 2) { return (v * 100).toFixed(dp) + '%'; }
function fmtNum(v)          { return Math.round(v).toLocaleString(); }

function fmtPval(p) {
  if (p < 0.0001) return '<0.0001';
  if (p > 0.9999) return '>0.9999';
  return p.toFixed(4);
}

function badge(label, type) {
  return `<span class="badge badge-${type}">${label}</span>`;
}

function alertBanner(type, icon, msg) {
  return `<div class="alert alert-${type}">
    <span class="alert-icon">${icon}</span>
    <span>${msg}</span>
  </div>`;
}

/* =====================================================================
   PART 1 — SAMPLE SIZE & DURATION
   ===================================================================== */

function calcPart1() {
  const cr    = parseFloat(el('p1-cr').value) / 100;
  const daily = parseFloat(el('p1-visitors').value);
  const mde   = parseFloat(el('p1-mde').value) / 100;
  const sig   = parseFloat(el('p1-sig').value);
  const power = parseFloat(el('p1-power').value);

  if (!cr || !daily || !mde || cr <= 0 || cr >= 1 || daily <= 0 || mde <= 0) return;

  const p1 = cr;
  const p2 = cr * (1 + mde);

  if (p2 >= 1) {
    setText('p1-n-per',    '—');
    setText('p1-n-total',  '—');
    setText('p1-days',     '—');
    setText('p1-weeks',    'MDE too large (target CR > 100%)');
    return;
  }

  const zAlpha = phiInv(sig);
  const zBeta  = power === 0.90 ? ZCRIT.p90 : ZCRIT.p80;

  // Two-proportion z-test sample size (one-sided)
  const n = Math.ceil(
    Math.pow(zAlpha + zBeta, 2) * (p1*(1-p1) + p2*(1-p2)) / Math.pow(p2 - p1, 2)
  );

  const dailyPerVar = daily / 2;
  const days  = Math.ceil(n / dailyPerVar);
  const weeks = (days / 7).toFixed(1);

  setText('p1-n-per',    n.toLocaleString());
  setText('p1-n-total',  (n * 2).toLocaleString());
  setText('p1-days',     days.toLocaleString());
  setText('p1-weeks',    `≈ ${weeks} weeks`);
  setText('p1-d-cr',     fmtPct(p1));
  setText('p1-d-target', fmtPct(p2));
  setText('p1-d-abs',    ((p2 - p1) * 100).toFixed(3) + ' pp');
  setText('p1-d-daily',  Math.floor(dailyPerVar).toLocaleString());
}

/* =====================================================================
   PART 2 — TEST EVALUATION
   ===================================================================== */

let testMode = 'standard';

function calcPart2() {
  const vA = parseFloat(el('p2-va').value);
  const cA = parseFloat(el('p2-ca').value);
  const vB = parseFloat(el('p2-vb').value);
  const cB = parseFloat(el('p2-cb').value);

  const alertsEl  = el('p2-alerts');
  const resultsEl = el('p2-results');
  const emptyEl   = el('p2-empty');

  alertsEl.innerHTML = '';

  const hasData = [vA, cA, vB, cB].every(v => !isNaN(v) && String(v) !== '');
  if (!hasData) {
    resultsEl.classList.add('hidden');
    emptyEl.classList.remove('hidden');
    return;
  }

  if (cA > vA || cB > vB) {
    alertsEl.innerHTML = alertBanner('error', '🚫', '<strong>Invalid data:</strong> Conversions cannot exceed visitors.');
    resultsEl.classList.add('hidden');
    emptyEl.classList.remove('hidden');
    return;
  }

  /* ── Core statistics ── */
  const crA    = cA / vA;
  const crB    = cB / vB;
  const pPool  = (cA + cB) / (vA + vB);
  const sePool = Math.sqrt(pPool * (1 - pPool) * (1/vA + 1/vB));
  const seA    = Math.sqrt(crA * (1 - crA) / vA);
  const seB    = Math.sqrt(crB * (1 - crB) / vB);
  const seDiff = Math.sqrt(seA*seA + seB*seB);

  const z      = sePool > 0 ? (crB - crA) / sePool : 0;
  const uplift = crA > 0 ? (crB - crA) / crA : 0;

  // P-values (one-sided)
  const pValStandard = 1 - phi(z);   // probability of improvement this large by chance
  const pValDnh      = phi(z);       // probability of harm this large by chance

  // Observed power against 95% threshold
  const obsPower = phi(Math.abs(z) - ZCRIT[95]) * 100;

  /* ── Alerts ── */
  let alerts = '';

  if (cA < 100 || cB < 100) {
    alerts += alertBanner('warn', '⚠️',
      '<strong>Low data warning:</strong> Fewer than 100 conversions in one or more variations. Results may be unreliable — consider collecting more data before drawing conclusions.');
  }

  // Sample Ratio Mismatch — chi-square test (1 df), p-value = 1 − erf(√(χ²/2))
  const total = vA + vB;
  const exp   = total / 2;
  const chiSq = (Math.pow(vA - exp, 2) + Math.pow(vB - exp, 2)) / exp;
  const srmP  = 1 - erf(Math.sqrt(chiSq / 2));
  if (srmP < 0.01) {
    alerts += alertBanner('error', '🚨',
      `<strong>Sample Ratio Mismatch (SRM) detected</strong> — visitor split is significantly unequal (p = ${srmP.toFixed(4)}). Your randomisation may be broken. Treat these results with caution until the root cause is investigated.`);
  }

  alertsEl.innerHTML = alerts;

  /* ── Summary metrics ── */
  setText('p2-cra', fmtPct(crA));
  setText('p2-crb', fmtPct(crB));

  const upliftEl = el('p2-uplift');
  upliftEl.textContent = (uplift >= 0 ? '+' : '') + (uplift * 100).toFixed(2) + '%';
  upliftEl.className   = 'metric-val ' + (uplift > 0 ? 'up' : uplift < 0 ? 'down' : '');

  setText('p2-z',     z.toFixed(4));
  setText('p2-power', Math.max(0, Math.min(100, obsPower)).toFixed(1) + '%');

  /* ── Mode-specific rendering ── */
  if (testMode === 'standard') {
    setText('p2-pval',     fmtPval(pValStandard));
    setText('p2-pval-lbl', 'P-Value');
    setText('p2-result-title', 'Standard A/B Test — Results by Confidence Level');
    setText('p2-result-desc',
      'Variant B is declared significant when the Z-Score exceeds the one-sided critical value. ' +
      'P-Value = probability of seeing this result (or stronger) by chance if there is no true difference.');
    setText('crit-lbl', 'Required Z (one-sided, improvement)');

    [95, 98, 99].forEach(lvl => {
      const zCrit = ZCRIT[lvl];
      const sig   = z > zCrit;
      setHTML(`r-${lvl}`,  badge(sig ? '✓ Significant' : '✗ Not Significant', sig ? 'success' : 'danger'));
      setText(`z-${lvl}`,  zCrit.toFixed(3));
      el(`yz-${lvl}`).textContent = z.toFixed(3);
      el(`yz-${lvl}`).className   = sig ? 'up' : '';
    });

  } else {
    // Do Not Harm
    // Harm is detected when B is significantly WORSE than A (z < −zCrit)
    setText('p2-pval',     fmtPval(pValDnh));
    setText('p2-pval-lbl', 'P-Value (harm)');
    setText('p2-result-title', 'Do Not Harm — Results by Confidence Level');
    setText('p2-result-desc',
      'Harm is detected when B is significantly worse than A (Z-Score falls below the lower critical threshold). ' +
      '"No Harm" means it is statistically safe to ship. P-Value (harm) = probability of B being at least this harmful by chance.');
    setText('crit-lbl', 'Lower critical Z (harm threshold)');

    [95, 98, 99].forEach(lvl => {
      const zCrit   = ZCRIT[lvl];
      const harmful = z < -zCrit;
      setHTML(`r-${lvl}`,  badge(harmful ? '✗ Harm Detected' : '✓ No Harm', harmful ? 'danger' : 'success'));
      setText(`z-${lvl}`,  (-zCrit).toFixed(3));
      el(`yz-${lvl}`).textContent = z.toFixed(3);
      el(`yz-${lvl}`).className   = harmful ? 'down' : 'up';
    });
  }

  /* ── Statistical details ── */
  setText('det-va',  vA.toLocaleString());
  setText('det-ca',  cA.toLocaleString());
  setText('det-cra', fmtPct(crA, 4));
  setText('det-vb',  vB.toLocaleString());
  setText('det-cb',  cB.toLocaleString());
  setText('det-crb', fmtPct(crB, 4));
  setText('det-sea', seA.toFixed(6));
  setText('det-seb', seB.toFixed(6));
  setText('det-sed', seDiff.toFixed(6));
  setText('det-z',   z.toFixed(6));

  resultsEl.classList.remove('hidden');
  emptyEl.classList.add('hidden');

  // Draw chart after results are visible (so canvas has layout width)
  requestAnimationFrame(() => drawChart(z, testMode));
}

// Redraw chart on window resize
let resizeTimer;
window.addEventListener('resize', () => {
  clearTimeout(resizeTimer);
  resizeTimer = setTimeout(() => {
    const canvas = el('p2-canvas');
    if (canvas && !el('p2-results').classList.contains('hidden')) {
      const z = parseFloat(el('p2-z').textContent);
      if (!isNaN(z)) drawChart(z, testMode);
    }
  }, 150);
});

/* =====================================================================
   DISTRIBUTION CHART
   ===================================================================== */

function drawChart(z, mode) {
  const canvas = el('p2-canvas');
  if (!canvas) return;

  const ctx = canvas.getContext('2d');
  const dpr = window.devicePixelRatio || 1;
  const cssW = canvas.parentElement.clientWidth - 48; // subtract card padding
  const cssH = 220;

  canvas.style.width  = cssW + 'px';
  canvas.style.height = cssH + 'px';
  canvas.width  = Math.round(cssW * dpr);
  canvas.height = Math.round(cssH * dpr);
  ctx.scale(dpr, dpr);

  const W = cssW, H = cssH;
  const padL = 24, padR = 24, padTop = 32, padBottom = 38;
  const plotW = W - padL - padR;
  const plotH = H - padTop - padBottom;
  const zMin = -4, zMax = 4;
  const maxPDF = normalPDF(0);

  function zToX(zv) {
    return padL + ((zv - zMin) / (zMax - zMin)) * plotW;
  }
  function pdfToY(pdf) {
    return padTop + plotH - (pdf / maxPDF) * plotH * 0.9;
  }

  // Fill area under curve between two z values
  function fillRegion(z1, z2, color) {
    const cz1 = Math.max(z1, zMin);
    const cz2 = Math.min(z2, zMax);
    if (cz1 >= cz2) return;
    const steps = 120;
    const dz = (cz2 - cz1) / steps;
    ctx.beginPath();
    ctx.moveTo(zToX(cz1), pdfToY(0));
    for (let i = 0; i <= steps; i++) {
      const zi = cz1 + i * dz;
      ctx.lineTo(zToX(zi), pdfToY(normalPDF(zi)));
    }
    ctx.lineTo(zToX(cz2), pdfToY(0));
    ctx.closePath();
    ctx.fillStyle = color;
    ctx.fill();
  }

  ctx.clearRect(0, 0, W, H);

  // Subtle vertical gridlines
  ctx.strokeStyle = '#e2e8f0';
  ctx.lineWidth = 1;
  [-3, -2, -1, 0, 1, 2, 3].forEach(v => {
    ctx.beginPath();
    ctx.moveTo(zToX(v), padTop);
    ctx.lineTo(zToX(v), padTop + plotH);
    ctx.stroke();
  });

  // Base bell curve fill
  fillRegion(zMin, zMax, 'rgba(219,234,254,0.35)');

  // Layered significance zones (lightest → darkest, outermost → innermost)
  if (mode === 'standard') {
    fillRegion(ZCRIT[95], zMax, 'rgba(134,239,172,0.50)');
    fillRegion(ZCRIT[98], zMax, 'rgba(74,222,128,0.55)');
    fillRegion(ZCRIT[99], zMax, 'rgba(22,163,74,0.60)');
  } else {
    fillRegion(zMin, -ZCRIT[95], 'rgba(252,165,165,0.50)');
    fillRegion(zMin, -ZCRIT[98], 'rgba(248,113,113,0.55)');
    fillRegion(zMin, -ZCRIT[99], 'rgba(220,38,38,0.60)');
  }

  // Bell curve outline
  ctx.beginPath();
  const steps = 300;
  for (let i = 0; i <= steps; i++) {
    const zi = zMin + (i / steps) * (zMax - zMin);
    if (i === 0) ctx.moveTo(zToX(zi), pdfToY(normalPDF(zi)));
    else         ctx.lineTo(zToX(zi), pdfToY(normalPDF(zi)));
  }
  ctx.strokeStyle = '#1d4ed8';
  ctx.lineWidth = 2;
  ctx.setLineDash([]);
  ctx.stroke();

  // X-axis baseline
  ctx.beginPath();
  ctx.moveTo(padL, padTop + plotH);
  ctx.lineTo(W - padR, padTop + plotH);
  ctx.strokeStyle = '#cbd5e1';
  ctx.lineWidth = 1;
  ctx.stroke();

  // Critical value dashed lines
  const critLines = mode === 'standard'
    ? [{ zv: ZCRIT[95], label: '95%' }, { zv: ZCRIT[98], label: '98%' }, { zv: ZCRIT[99], label: '99%' }]
    : [{ zv: -ZCRIT[99], label: '99%' }, { zv: -ZCRIT[98], label: '98%' }, { zv: -ZCRIT[95], label: '95%' }];

  const critColors = ['#9ca3af', '#6b7280', '#374151'];
  ctx.font = `10px -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif`;

  critLines.forEach(({ zv, label }, i) => {
    const x = zToX(zv);
    ctx.setLineDash([4, 3]);
    ctx.beginPath();
    ctx.moveTo(x, padTop);
    ctx.lineTo(x, padTop + plotH);
    ctx.strokeStyle = critColors[i];
    ctx.lineWidth = 1.5;
    ctx.stroke();
    ctx.setLineDash([]);

    ctx.fillStyle = critColors[i];
    ctx.textAlign = 'center';
    ctx.fillText(label, x, padTop - 10);
  });

  // Observed z-score vertical line
  const clampedZ = Math.min(Math.max(z, zMin + 0.01), zMax - 0.01);
  const zX = zToX(clampedZ);

  const lineColor = mode === 'standard'
    ? (z > ZCRIT[95] ? '#15803d' : '#1d4ed8')
    : (z < -ZCRIT[95] ? '#b91c1c' : '#15803d');

  ctx.beginPath();
  ctx.moveTo(zX, padTop - 4);
  ctx.lineTo(zX, padTop + plotH);
  ctx.strokeStyle = lineColor;
  ctx.lineWidth = 2.5;
  ctx.setLineDash([]);
  ctx.stroke();

  // Arrow head at top of z-score line
  ctx.beginPath();
  ctx.moveTo(zX, padTop - 4);
  ctx.lineTo(zX - 4, padTop + 6);
  ctx.lineTo(zX + 4, padTop + 6);
  ctx.closePath();
  ctx.fillStyle = lineColor;
  ctx.fill();

  // Z-score label at bottom
  ctx.font = `bold 11px -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif`;
  ctx.fillStyle = lineColor;
  ctx.textAlign = 'center';
  // Offset label slightly if it overlaps a critical value line
  const labelY = padTop + plotH + 14;
  ctx.fillText(`z = ${z.toFixed(3)}`, zX, labelY);

  // X-axis tick labels
  ctx.font = `10px -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif`;
  ctx.fillStyle = '#94a3b8';
  [-4, -3, -2, -1, 0, 1, 2, 3, 4].forEach(v => {
    ctx.textAlign = 'center';
    ctx.fillText(v, zToX(v), padTop + plotH + 26);
  });

  // Update legend HTML
  if (mode === 'standard') {
    setText('p2-chart-desc', 'Green shaded regions show where Variant B would be statistically significant. Your Z-score (arrow) needs to land in a shaded zone to reach that confidence level.');
    setHTML('p2-chart-legend', `
      <div class="legend-item"><span class="legend-swatch" style="background:rgba(134,239,172,0.7)"></span>Significant at 95%</div>
      <div class="legend-item"><span class="legend-swatch" style="background:rgba(74,222,128,0.8)"></span>Significant at 98%</div>
      <div class="legend-item"><span class="legend-swatch" style="background:rgba(22,163,74,0.85)"></span>Significant at 99%</div>
      <div class="legend-item"><span class="legend-swatch" style="background:${lineColor};opacity:.8"></span>Your Z-score (${z.toFixed(3)})</div>
    `);
  } else {
    setText('p2-chart-desc', 'Red shaded regions show where Variant B would be causing statistically significant harm. Your Z-score (arrow) must stay outside the red zones to pass Do Not Harm.');
    setHTML('p2-chart-legend', `
      <div class="legend-item"><span class="legend-swatch" style="background:rgba(252,165,165,0.7)"></span>Harm detected at 95%</div>
      <div class="legend-item"><span class="legend-swatch" style="background:rgba(248,113,113,0.8)"></span>Harm detected at 98%</div>
      <div class="legend-item"><span class="legend-swatch" style="background:rgba(220,38,38,0.85)"></span>Harm detected at 99%</div>
      <div class="legend-item"><span class="legend-swatch" style="background:${lineColor};opacity:.8"></span>Your Z-score (${z.toFixed(3)})</div>
    `);
  }
}

/* =====================================================================
   EVENT LISTENERS
   ===================================================================== */

// Tab switching
document.querySelectorAll('.tab-btn').forEach(btn => {
  btn.addEventListener('click', () => {
    document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
    document.querySelectorAll('.tab-pane').forEach(p => p.classList.remove('active'));
    btn.classList.add('active');
    el(btn.dataset.tab).classList.add('active');
  });
});

// Part 1 — live calculation
['p1-cr', 'p1-visitors', 'p1-mde', 'p1-sig', 'p1-power'].forEach(id => {
  el(id).addEventListener('input', calcPart1);
});

// Part 2 — live calculation
['p2-va', 'p2-ca', 'p2-vb', 'p2-cb'].forEach(id => {
  el(id).addEventListener('input', calcPart2);
});

// Test type toggle
document.querySelectorAll('.toggle-btn').forEach(btn => {
  btn.addEventListener('click', () => {
    document.querySelectorAll('.toggle-btn').forEach(b => b.classList.remove('active'));
    btn.classList.add('active');
    testMode = btn.dataset.mode;

    el('p2-mode-hint').textContent = testMode === 'standard'
      ? 'Tests whether Variant B shows a statistically significant improvement over Control A.'
      : 'Tests whether Variant B causes no significant harm to conversion rates. Use when shipping a change without needing to prove uplift — you just want to confirm it is not hurting performance.';

    calcPart2();
  });
});

// Initial calculation
calcPart1();
