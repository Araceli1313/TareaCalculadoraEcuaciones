document.addEventListener('DOMContentLoaded', () => {
    const eqInput = document.getElementById('equationInput');
    const x0Input = document.getElementById('x0Input');
    const y0Input = document.getElementById('y0Input');
    const xEndInput = document.getElementById('xEndInput');
    const hInput = document.getElementById('hInput');
    const calcBtn = document.getElementById('calculateBtn');
    const solutionText = document.getElementById('solutionText');
    const errorText = document.getElementById('errorText');
    const chartCanvas = document.getElementById('solutionChart');
    const downloadBtn = document.getElementById('downloadChart');
    const methodSelect = document.getElementById('methodSelect');

    let chartInstance = null;

    function parseRHS(equation) {
      if (!equation) return null;
      const idx = equation.indexOf('=');
      let rhs = idx >= 0 ? equation.slice(idx + 1) : equation;
      rhs = rhs.replace(/y\'/g, '');
      rhs = rhs.replace(/dy\/dx/gi, '');
      rhs = rhs.replace(/\^/g, '**');
      return rhs.trim();
    }

    function makeFunction(expr) {
      try {
        const node = math.parse(expr);
        const compiled = node.compile();
        return (x, y) => {
          const scope = { x: Number(x), y: Number(y) };
          const val = compiled.evaluate(scope);
          return Number(val);
        };
      } catch (err) {
        throw new Error('No se pudo parsear la expresión: ' + err.message);
      }
    }

    function solveRK4(f, x0, y0, xEnd, h) {
      const xs = [x0];
      const ys = [y0];
      let x = x0;
      let y = y0;
      const n = Math.max(1, Math.ceil(Math.abs((xEnd - x0) / h)));
      for (let i = 0; i < n; i++) {
        const hstep = (xEnd >= x0) ? Math.min(h, xEnd - x) : Math.max(-h, xEnd - x);
        const k1 = f(x, y);
        const k2 = f(x + hstep / 2, y + (hstep * k1) / 2);
        const k3 = f(x + hstep / 2, y + (hstep * k2) / 2);
        const k4 = f(x + hstep, y + hstep * k3);
        y = y + (hstep / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
        x = x + hstep;
        xs.push(x);
        ys.push(y);
        if ((hstep > 0 && x >= xEnd) || (hstep < 0 && x <= xEnd)) break;
      }
      return { xs, ys };
    }

    function solveLinearByIntegratingFactor(exprStr, x0, y0, xEnd, h) {
      const node = math.parse(exprStr);
      const dY = math.derivative(node, 'y');
      const dYStr = dY.toString();
      if (dYStr.includes('y')) {
        throw new Error('La ecuación no parece lineal (el coeficiente depende de y).');
      }
      const aComp = dY.compile();
      const fCompiled = node.compile();

      const aFunc = x => {
        return Number(aComp.evaluate({ x, y: 0 }));
      };
      const bFunc = x => {
        return Number(fCompiled.evaluate({ x, y: 0 }));
      };

      const xs = [x0];
      const ys = [y0];
      let x = x0;
      let n = Math.max(1, Math.ceil(Math.abs((xEnd - x0) / h)));
      let cumA = 0;
      let cumMuB = 0;
      const mu0 = Math.exp(-0);
      let prevA = aFunc(x0);
      let prevMuB = Math.exp(-cumA) * bFunc(x0);

      for (let i = 0; i < n; i++) {
        const hstep = (xEnd >= x0) ? Math.min(h, xEnd - x) : Math.max(-h, xEnd - x);
        const xNext = x + hstep;
        const aNext = aFunc(xNext);
        cumA += 0.5 * (prevA + aNext) * hstep;
        const mu = Math.exp(-cumA);
        const muPrev = Math.exp(- (cumA - 0.5 * (prevA + aNext) * hstep));
        const bNext = bFunc(xNext);
        const muBNext = mu * bNext;
        cumMuB += 0.5 * (prevMuB + muBNext) * hstep;

        const y = (1 / mu) * (mu0 * y0 + cumMuB);
        xs.push(xNext);
        ys.push(y);

        x = xNext;
        prevA = aNext;
        prevMuB = muBNext;
        if ((hstep > 0 && x >= xEnd) || (hstep < 0 && x <= xEnd)) break;
      }
      return { xs, ys };
    }

    function solveHomogeneous(exprStr, x0, y0, xEnd, h) {
      const node = math.parse(exprStr);
      const compiled = node.compile();
      const testV = 2.5;
      const val1 = compiled.evaluate({ x: 1.234, y: testV * 1.234 });
      const val2 = compiled.evaluate({ x: 2.345, y: testV * 2.345 });
      if (Math.abs(val1 - val2) > 1e-8) {
        throw new Error('No parece ser homogénea (no es reducible por y = v*x).');
      }

      const fVnode = math.parse(exprStr.replace(/y/g, '(v*x)'));
      const fVcompiled = fVnode.compile();

      function dv_dx(x, v) {
        const Fv = Number(fVcompiled.evaluate({ x, v }));
        return (Fv - v) / x;
      }

      const xs = [x0];
      const ys = [y0];
      let x = x0;
      let v = (x0 === 0) ? (y0 / (1e-8)) : (y0 / x0);
      const n = Math.max(1, Math.ceil(Math.abs((xEnd - x0) / h)));
      for (let i = 0; i < n; i++) {
        const hstep = (xEnd >= x0) ? Math.min(h, xEnd - x) : Math.max(-h, xEnd - x);
        const k1 = dv_dx(x === 0 ? (x + 1e-8) : x, v);
        const k2 = dv_dx(x + hstep / 2, v + (hstep * k1) / 2);
        const k3 = dv_dx(x + hstep / 2, v + (hstep * k2) / 2);
        const k4 = dv_dx(x + hstep, v + hstep * k3);
        v = v + (hstep / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
        x = x + hstep;
        const y = v * x;
        xs.push(x);
        ys.push(y);
        if ((hstep > 0 && x >= xEnd) || (hstep < 0 && x <= xEnd)) break;
      }
      return { xs, ys };
    }

    function solveBernoulli(exprStr, x0, y0, xEnd, h) {
      const node = math.parse(exprStr);
      const compiled = node.compile();

      let n = null;
      node.traverse(function (node2, path, parent) {
        if (node2.type === 'OperatorNode' && node2.op === '^') {
          const [base, exp] = node2.args;
          if (base && base.isSymbolNode && base.name === 'y' && exp.isConstantNode) {
            n = Number(exp.value);
          }
        }
      });
      if (n === null || n === 1) {
        throw new Error('No se detectó un término tipo y^n apropiado para Bernoulli.');
      }

      const dY = math.derivative(node, 'y').compile();
      const aFunc = x => Number(dY.evaluate({ x, y: 0 }));
      const fCompiled = compiled;
      const bFunc = x => Number(fCompiled.evaluate({ x, y: 1 }) - aFunc(x));

      const alpha = x => (1 - n) * aFunc(x);
      const beta = x => (1 - n) * bFunc(x);

      const xs = [x0];
      const zs = [Math.pow(y0, 1 - n)];
      let x = x0;
      let cumAlpha = 0;
      let cumMuBeta = 0;
      let prevAlpha = alpha(x0);
      let prevMuBeta = Math.exp(-cumAlpha) * beta(x0);
      const mu0 = Math.exp(-0);
      const nSteps = Math.max(1, Math.ceil(Math.abs((xEnd - x0) / h)));
      for (let i = 0; i < nSteps; i++) {
        const hstep = (xEnd >= x0) ? Math.min(h, xEnd - x) : Math.max(-h, xEnd - x);
        const xNext = x + hstep;
        const alphaNext = alpha(xNext);
        cumAlpha += 0.5 * (prevAlpha + alphaNext) * hstep;
        const mu = Math.exp(-cumAlpha);
        const betaNext = beta(xNext);
        const muBetaNext = mu * betaNext;
        cumMuBeta += 0.5 * (prevMuBeta + muBetaNext) * hstep;
        const z = (1 / mu) * (mu0 * zs[0] + cumMuBeta);
        const y = Math.pow(z, 1 / (1 - n));
        xs.push(xNext);
        zs.push(z);
        x = xNext;
        prevAlpha = alphaNext;
        prevMuBeta = muBetaNext;
        if ((hstep > 0 && x >= xEnd) || (hstep < 0 && x <= xEnd)) break;
      }
      const ysRes = zs.map(z => Math.pow(z, 1 / (1 - n)));
      return { xs, ys: ysRes };
    }

    function isRiccati(node) {
      let maxPow = 0;
      node.traverse(function (n) {
        if (n.type === 'OperatorNode' && n.op === '^') {
          const [base, exp] = n.args;
          if (base.isSymbolNode && base.name === 'y' && exp.isConstantNode) {
            maxPow = Math.max(maxPow, Number(exp.value));
          }
        }
      });
      return maxPow >= 2;
    }

    function renderChart(xs, ys) {
      if (chartInstance) {
        chartInstance.destroy();
        chartInstance = null;
      }
      const labels = xs.map(v => Number(v.toFixed(6)));
      const data = {
        labels,
        datasets: [
          {
            label: 'y(x)',
            data: ys,
            borderColor: '#2563eb',
            backgroundColor: 'rgba(37,99,235,0.08)',
            pointRadius: 0,
            borderWidth: 2,
            tension: 0.2,
          },
        ],
      };
      chartInstance = new Chart(chartCanvas.getContext('2d'), {
        type: 'line',
        data,
        options: {
          responsive: true,
          maintainAspectRatio: false,
          scales: {
            x: { display: true, title: { display: true, text: 'x' } },
            y: { display: true, title: { display: true, text: 'y' } },
          },
        },
      });
    }

    calcBtn.addEventListener('click', () => {
      errorText.textContent = '';
      try {
        const eq = eqInput.value.trim();
        if (!eq) throw new Error('Introduce una ecuación, por ejemplo: y\' = 2*x + y');
        const rhs = parseRHS(eq);
        if (!rhs) throw new Error('No se encontró el lado derecho de la ecuación.');
        const method = (methodSelect && methodSelect.value) || 'linear';
        const x0 = Number(x0Input.value);
        const y0 = Number(y0Input.value);
        const xEnd = Number(xEndInput.value);
        let h = Number(hInput.value);
        if (!isFinite(h) || h === 0) h = 0.01;
        if (xEnd < x0 && h > 0) h = -h;

        let result = null;
        try {
          if (method === 'linear') {
            result = solveLinearByIntegratingFactor(rhs, x0, y0, xEnd, h);
          } else if (method === 'homogeneous') {
            result = solveHomogeneous(rhs, x0, y0, xEnd, h);
          } else if (method === 'bernoulli') {
            result = solveBernoulli(rhs, x0, y0, xEnd, h);
          } else if (method === 'riccati') {
            const node = math.parse(rhs);
            if (!isRiccati(node)) {
              throw new Error('No se detectó una forma tipo Riccati (no se encontró término y^2).');
            }
            const f = makeFunction(rhs);
            result = solveRK4(f, x0, y0, xEnd, h);
          } else {
            const f = makeFunction(rhs);
            result = solveRK4(f, x0, y0, xEnd, h);
          }
        } catch (errMethod) {
          console.warn('Método especializado falló, usando RK4:', errMethod.message);
          const f = makeFunction(rhs);
          result = solveRK4(f, x0, y0, xEnd, h);
          errorText.textContent = 'Método especializado falló: ' + errMethod.message + ' — se usó RK4 como fallback.';
        }

        const { xs, ys } = result;
        function escapeHtml(unsafe) {
          return unsafe
            .replace(/&/g, '&amp;')
            .replace(/</g, '&lt;')
            .replace(/>/g, '&gt;')
            .replace(/"/g, '&quot;')
            .replace(/'/g, '&#039;');
        }

        const methodName = (methodSelect && methodSelect.selectedOptions && methodSelect.selectedOptions[0])
          ? methodSelect.selectedOptions[0].textContent
          : (method || 'Numérico (RK4)');

        const eqEscaped = escapeHtml(eq);
        const finalVal = ys[ys.length - 1];
        const finalX = xs[xs.length - 1];
        solutionText.innerHTML = `
          <div class="flex flex-col gap-2">
            <div class="text-sm text-gray-600 dark:text-gray-400">Método: <span class="font-medium">${escapeHtml(methodName)}</span></div>
            <div class="font-mono text-sm text-gray-800 dark:text-gray-200">Ecuación: <code>${eqEscaped}</code></div>
            <div class="text-base font-semibold mt-2">Valor final: y(${Number(finalX).toPrecision(6)}) ≈ ${Number(finalVal).toPrecision(6)}</div>
          </div>
        `;

        renderChart(xs, ys);
        try {
          solutionText.scrollIntoView({ behavior: 'smooth', block: 'center' });
        } catch (e) {
        }
      } catch (err) {
        errorText.textContent = err.message || String(err);
      }
    });

    downloadBtn.addEventListener('click', () => {
      if (!chartInstance) return;
      const url = chartInstance.toBase64Image();
      const a = document.createElement('a');
      a.href = url;
      a.download = 'solucion.png';
      a.click();
    });

    const copyBtn = document.querySelector('button[aria-label="Copiar Solución"]');
    if (copyBtn) {
      copyBtn.addEventListener('click', async () => {
        try {
          await navigator.clipboard.writeText(solutionText.textContent || '');
          copyBtn.textContent = 'Copiado';
          setTimeout(() => (copyBtn.innerHTML = '<span class="material-symbols-outlined text-base">content_copy</span>'), 1200);
        } catch (e) {
        }
      });
    }
  });