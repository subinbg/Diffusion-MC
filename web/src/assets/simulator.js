// DMC Web Interface - Main Application

import init, { WebSimulation } from './pkg/web.js';

// Global state
let simulation = null;
let animationId = null;
let isRunning = false;
let energyChart = null;
let populationChart = null;
let currentColorTheme = 'magma';
let lastMaxDensity = 0;

// Known bond lengths (in Bohr radii)
const BOND_LENGTHS = {
    hydrogen: 0,       // Single atom
    h2_ion: 2.0,       // H2+ equilibrium bond length
    h2_molecule: 1.4   // H2 equilibrium bond length
};

// Exact energies (in Hartree)
const EXACT_ENERGIES = {
    hydrogen: -0.5,
    h2_ion: -0.6026,
    h2_molecule: -1.1745
};

// Color map definitions
const COLOR_MAPS = {
    viridis: [
        [68, 1, 84], [72, 35, 116], [64, 67, 135], [52, 94, 141],
        [41, 120, 142], [32, 144, 140], [34, 167, 132], [68, 190, 112],
        [121, 209, 81], [189, 222, 38], [253, 231, 37]
    ],
    plasma: [
        [13, 8, 135], [75, 3, 161], [125, 3, 168], [168, 34, 150],
        [203, 70, 121], [229, 107, 93], [248, 148, 65], [253, 195, 40],
        [240, 249, 33]
    ],
    inferno: [
        [0, 0, 4], [22, 11, 57], [66, 10, 104], [106, 23, 110],
        [147, 38, 103], [188, 55, 84], [221, 81, 58], [243, 118, 27],
        [252, 166, 10], [252, 220, 38], [252, 255, 164]
    ],
    magma: [
        [0, 0, 4], [20, 14, 54], [57, 15, 110], [102, 17, 136],
        [143, 36, 143], [183, 55, 137], [218, 82, 123], [241, 120, 108],
        [253, 166, 122], [254, 214, 164], [252, 253, 191]
    ],
    grayscale: [
        [0, 0, 0], [28, 28, 28], [57, 57, 57], [85, 85, 85],
        [113, 113, 113], [142, 142, 142], [170, 170, 170], [198, 198, 198],
        [227, 227, 227], [255, 255, 255]
    ]
};

// DOM elements
const elements = {
    system: document.getElementById('system'),
    algorithm: document.getElementById('algorithm'),
    numWalkers: document.getElementById('num_walkers'),
    timeStep: document.getElementById('time_step'),
    totalSteps: document.getElementById('total_steps'),
    equilibrationSteps: document.getElementById('equilibration_steps'),
    updateInterval: document.getElementById('update_interval'),
    updateIntervalValue: document.getElementById('update_interval_value'),
    startBtn: document.getElementById('start_btn'),
    pauseBtn: document.getElementById('pause_btn'),
    resetBtn: document.getElementById('reset_btn'),
    progress: document.getElementById('progress'),
    statStep: document.getElementById('stat_step'),
    statTotal: document.getElementById('stat_total'),
    statPhase: document.getElementById('stat_phase'),
    statEnergy: document.getElementById('stat_energy'),
    statError: document.getElementById('stat_error'),
    statExact: document.getElementById('stat_exact'),
    statPopulation: document.getElementById('stat_population'),
    statAcceptance: document.getElementById('stat_acceptance'),
    statDeviation: document.getElementById('stat_deviation'),
    densityPlane: document.getElementById('density_plane'),
    slicePosition: document.getElementById('slice_position'),
    sliceValue: document.getElementById('slice_value'),
    densityCanvas: document.getElementById('density_canvas'),
    colorTheme: document.getElementById('color_theme'),
    colorbar: document.getElementById('colorbar'),
    colorbarMax: document.getElementById('colorbar_max'),
    colorbarMin: document.getElementById('colorbar_min'),
    densityInfo: document.getElementById('density_info'),
};

// Initialize application
async function main() {
    await init();
    setupEventListeners();
    initCharts();
    updateUIState();
    updateColorbar();
    console.log('DMC Web Interface initialized');
}

// Set up event listeners
function setupEventListeners() {
    elements.system.addEventListener('change', onSystemChange);
    elements.updateInterval.addEventListener('input', onUpdateIntervalChange);
    elements.startBtn.addEventListener('click', onStart);
    elements.pauseBtn.addEventListener('click', onPause);
    elements.resetBtn.addEventListener('click', onReset);
    elements.densityPlane.addEventListener('change', updateDensity);
    elements.slicePosition.addEventListener('input', onSliceChange);
    elements.colorTheme.addEventListener('change', onColorThemeChange);
}

// System change handler
function onSystemChange() {
    const system = elements.system.value;
    elements.statExact.textContent = EXACT_ENERGIES[system].toFixed(4);
}

// Update interval change handler
function onUpdateIntervalChange() {
    const value = elements.updateInterval.value;
    elements.updateIntervalValue.textContent = value;
    if (simulation) {
        simulation.set_update_interval(parseInt(value));
    }
}

// Slice position change handler
function onSliceChange() {
    elements.sliceValue.textContent = parseFloat(elements.slicePosition.value).toFixed(1);
    if (simulation && !simulation.is_finished()) {
        updateDensity();
    }
}

// Color theme change handler
function onColorThemeChange() {
    currentColorTheme = elements.colorTheme.value;
    updateColorbar();
    if (simulation) {
        updateDensity();
    }
}

// Start simulation
function onStart() {
    if (!simulation) {
        createSimulation();
    }

    if (simulation) {
        isRunning = true;
        updateUIState();
        runSimulation();
    }
}

// Pause simulation
function onPause() {
    isRunning = false;
    if (animationId) {
        cancelAnimationFrame(animationId);
        animationId = null;
    }
    updateUIState();
}

// Reset simulation
function onReset() {
    onPause();
    simulation = null;
    resetCharts();
    resetStats();
    clearDensity();
    updateUIState();
}

// Create new simulation
function createSimulation() {
    const system = elements.system.value;
    const config = {
        system: system,
        algorithm: elements.algorithm.value,
        num_walkers: parseInt(elements.numWalkers.value),
        time_step: parseFloat(elements.timeStep.value),
        total_steps: parseInt(elements.totalSteps.value),
        equilibration_steps: parseInt(elements.equilibrationSteps.value),
        seed: Math.floor(Math.random() * 1000000),
    };

    // Use fixed bond lengths
    if (system !== 'hydrogen') {
        config.bond_length = BOND_LENGTHS[system];
    }

    if (elements.algorithm.value === 'importance_sampled') {
        config.alpha = 1.0;
        if (system === 'h2_molecule') {
            config.alpha = 1.2;
            config.jastrow_b = 0.5;
        }
    }

    try {
        simulation = new WebSimulation(JSON.stringify(config));
        simulation.set_update_interval(parseInt(elements.updateInterval.value));
        elements.statTotal.textContent = config.total_steps;
        console.log('Simulation created:', simulation.system_name());
    } catch (e) {
        console.error('Failed to create simulation:', e);
        alert('Failed to create simulation: ' + e);
    }
}

// Run simulation loop
function runSimulation() {
    if (!isRunning || !simulation || simulation.is_finished()) {
        if (simulation && simulation.is_finished()) {
            isRunning = false;
            updateUIState();
            updateStats();
            updateDensity();
        }
        return;
    }

    const interval = parseInt(elements.updateInterval.value);

    try {
        const result = simulation.run_steps(interval);

        if (result) {
            updateStats(result);
            updateCharts();

            // Update density less frequently
            if (simulation.current_step() % (interval * 10) === 0) {
                updateDensity();
            }
        }

        animationId = requestAnimationFrame(runSimulation);
    } catch (e) {
        console.error('Simulation error:', e);
        onPause();
    }
}

// Update UI state
function updateUIState() {
    const hasSimulation = simulation !== null;
    const finished = simulation ? simulation.is_finished() : false;

    elements.startBtn.disabled = isRunning || finished;
    elements.pauseBtn.disabled = !isRunning;
    elements.resetBtn.disabled = !hasSimulation && !isRunning;

    // Disable config inputs while running
    const configDisabled = isRunning;
    elements.system.disabled = configDisabled;
    elements.algorithm.disabled = configDisabled;
    elements.numWalkers.disabled = configDisabled;
    elements.timeStep.disabled = configDisabled;
    elements.totalSteps.disabled = configDisabled;
    elements.equilibrationSteps.disabled = configDisabled;
}

// Update statistics display
function updateStats(stepResult) {
    if (!simulation) return;

    const step = simulation.current_step();
    const total = simulation.total_steps();
    const isEquilibrated = simulation.is_equilibrated();
    const finished = simulation.is_finished();

    elements.statStep.textContent = step;
    elements.progress.value = (step / total) * 100;

    // Phase
    let phase;
    if (finished) {
        phase = 'Complete';
        elements.statPhase.className = 'stat-value phase-complete';
    } else if (isEquilibrated) {
        phase = 'Sampling';
        elements.statPhase.className = 'stat-value phase-sampling';
    } else {
        phase = 'Equilibrating';
        elements.statPhase.className = 'stat-value phase-equilibrating';
    }
    elements.statPhase.textContent = phase;

    // Population
    elements.statPopulation.textContent = simulation.population_size();

    // Acceptance ratio:
    // - Pure DMC: no Metropolis step, so acceptance_ratio is always null -> "N/A"
    // - Importance-sampled: has acceptance_ratio -> show percentage
    if (elements.algorithm.value === 'pure') {
        elements.statAcceptance.textContent = 'N/A';
    } else if (stepResult && stepResult.acceptance_ratio != null && !isNaN(stepResult.acceptance_ratio)) {
        elements.statAcceptance.textContent = (stepResult.acceptance_ratio * 100).toFixed(1);
    }

    // Get statistics
    try {
        const stats = simulation.get_statistics();
        if (stats && stats.sample_count > 0) {
            elements.statEnergy.textContent = stats.mean_energy.toFixed(6);
            elements.statError.textContent = stats.energy_error.toFixed(6);

            // Deviation from exact
            if (stats.exact_energy != null) {
                const deviation = Math.abs((stats.mean_energy - stats.exact_energy) / stats.exact_energy * 100);
                elements.statDeviation.textContent = deviation.toFixed(2);

                // Color code deviation
                let devClass = 'deviation-good';
                if (deviation > 2) devClass = 'deviation-bad';
                else if (deviation > 0.5) devClass = 'deviation-warning';
                elements.statDeviation.className = 'stat-value ' + devClass;
            }
        }
    } catch (e) {
        console.warn('Stats error:', e);
    }
}

// Reset statistics display
function resetStats() {
    elements.statStep.textContent = '0';
    elements.statPhase.textContent = 'Ready';
    elements.statPhase.className = 'stat-value phase-ready';
    elements.statEnergy.textContent = '-';
    elements.statError.textContent = '-';
    elements.statPopulation.textContent = '0';
    elements.statAcceptance.textContent = '-';
    elements.statDeviation.textContent = '-';
    elements.statDeviation.className = 'stat-value';
    elements.progress.value = 0;
    elements.densityInfo.textContent = '';
    lastMaxDensity = 0;
    updateColorbar();
}

// Initialize charts
function initCharts() {
    const chartOptions = {
        responsive: true,
        maintainAspectRatio: true,
        animation: false,
        scales: {
            x: {
                title: { display: true, text: 'Step' }
            }
        },
        plugins: {
            legend: { display: false }
        }
    };

    energyChart = new Chart(document.getElementById('energy_chart'), {
        type: 'line',
        data: {
            labels: [],
            datasets: [{
                label: 'Energy',
                data: [],
                borderColor: '#0077b6',
                borderWidth: 1,
                pointRadius: 0,
                tension: 0.1
            }, {
                label: 'Reference',
                data: [],
                borderColor: '#d62828',
                borderWidth: 2,
                borderDash: [5, 5],
                pointRadius: 0
            }]
        },
        options: {
            ...chartOptions,
            scales: {
                ...chartOptions.scales,
                y: { title: { display: true, text: 'Energy (Ha)' } }
            },
            plugins: {
                legend: { display: true, position: 'top' }
            }
        }
    });

    populationChart = new Chart(document.getElementById('population_chart'), {
        type: 'line',
        data: {
            labels: [],
            datasets: [{
                label: 'Population',
                data: [],
                borderColor: '#2a9d8f',
                borderWidth: 1,
                pointRadius: 0,
                fill: true,
                backgroundColor: 'rgba(42, 157, 143, 0.1)'
            }]
        },
        options: {
            ...chartOptions,
            scales: {
                ...chartOptions.scales,
                y: { title: { display: true, text: 'Walkers' } }
            }
        }
    });
}

// Update charts
function updateCharts() {
    if (!simulation) return;

    try {
        const energyTrace = simulation.get_energy_trace();
        const populationTrace = simulation.get_population_trace();

        if (energyTrace.length === 0) return;

        // Downsample for performance
        const maxPoints = 500;
        const step = Math.max(1, Math.floor(energyTrace.length / maxPoints));

        const labels = [];
        const energyData = [];
        const populationData = [];
        const referenceData = [];

        const exactEnergy = simulation.exact_energy();
        const equilibration = simulation.equilibration_steps();

        for (let i = 0; i < energyTrace.length; i += step) {
            labels.push(equilibration + i);
            energyData.push(energyTrace[i]);
            populationData.push(populationTrace[i]);
            if (exactEnergy != null) {
                referenceData.push(exactEnergy);
            }
        }

        energyChart.data.labels = labels;
        energyChart.data.datasets[0].data = energyData;
        energyChart.data.datasets[1].data = referenceData;
        energyChart.update('none');

        populationChart.data.labels = labels;
        populationChart.data.datasets[0].data = populationData;
        populationChart.update('none');
    } catch (e) {
        console.warn('Chart update error:', e);
    }
}

// Reset charts
function resetCharts() {
    energyChart.data.labels = [];
    energyChart.data.datasets[0].data = [];
    energyChart.data.datasets[1].data = [];
    energyChart.update();

    populationChart.data.labels = [];
    populationChart.data.datasets[0].data = [];
    populationChart.update();
}

// Get color from colormap
function getColor(t, colormap) {
    t = Math.max(0, Math.min(1, t));
    const colors = COLOR_MAPS[colormap] || COLOR_MAPS.viridis;
    const n = colors.length - 1;
    const idx = t * n;
    const i = Math.floor(idx);
    const f = idx - i;

    if (i >= n) return colors[n];
    if (i < 0) return colors[0];

    const c1 = colors[i];
    const c2 = colors[i + 1];

    return [
        Math.round(c1[0] + f * (c2[0] - c1[0])),
        Math.round(c1[1] + f * (c2[1] - c1[1])),
        Math.round(c1[2] + f * (c2[2] - c1[2]))
    ];
}

// Update colorbar gradient
function updateColorbar() {
    // Generate smooth gradient using the same interpolation as the density figure
    // CSS gradient stops must be in ascending order (0% at top, 100% at bottom)
    const numStops = 50;
    const stops = [];

    for (let i = 0; i <= numStops; i++) {
        const pct = (i / numStops) * 100;  // 0% to 100% (top to bottom)
        const t = 1 - (i / numStops);       // 1 to 0 (high density at top, low at bottom)
        const color = getColor(t, currentColorTheme);
        stops.push(`rgb(${color[0]}, ${color[1]}, ${color[2]}) ${pct.toFixed(1)}%`);
    }

    elements.colorbar.style.background = `linear-gradient(to bottom, ${stops.join(', ')})`;

    // Update colorbar labels
    if (lastMaxDensity > 0) {
        elements.colorbarMax.textContent = lastMaxDensity.toExponential(2);
    } else {
        elements.colorbarMax.textContent = '1.0';
    }
    elements.colorbarMin.textContent = '0.0';
}

// Update density visualization
function updateDensity() {
    if (!simulation) return;

    const canvas = elements.densityCanvas;
    const ctx = canvas.getContext('2d');
    const size = canvas.width;

    try {
        const plane = elements.densityPlane.value;
        const slicePos = parseFloat(elements.slicePosition.value);
        const gridSize = 100;
        const bounds = 5.0;
        const sigma = 0.3;

        const density = simulation.get_density_2d(plane, slicePos, gridSize, bounds, sigma);

        // Find max density for normalization
        let maxDensity = 0;
        for (const d of density) {
            if (d > maxDensity) maxDensity = d;
        }

        if (maxDensity === 0) {
            clearDensity();
            return;
        }

        lastMaxDensity = maxDensity;
        updateColorbar();

        // Create image data
        const imageData = ctx.createImageData(size, size);
        const scaleX = size / gridSize;
        const scaleY = size / gridSize;

        for (let y = 0; y < size; y++) {
            for (let x = 0; x < size; x++) {
                const gridX = Math.floor(x / scaleX);
                const gridY = Math.floor(y / scaleY);
                const idx = gridY * gridSize + gridX;
                const value = density[idx] / maxDensity;

                const color = getColor(value, currentColorTheme);
                const pixelIdx = (y * size + x) * 4;
                imageData.data[pixelIdx] = color[0];
                imageData.data[pixelIdx + 1] = color[1];
                imageData.data[pixelIdx + 2] = color[2];
                imageData.data[pixelIdx + 3] = 255;
            }
        }

        ctx.putImageData(imageData, 0, 0);

        // Draw nuclear positions and bond length
        drawNucleiAndBond(ctx, plane, slicePos, bounds, size);

        // Update density info
        updateDensityInfo(plane, slicePos);

    } catch (e) {
        console.warn('Density update error:', e);
    }
}

// Update density info text
function updateDensityInfo(plane, slicePos) {
    const system = elements.system.value;
    const bondLength = BOND_LENGTHS[system];

    let info = `Plane: ${plane.toUpperCase()}, Slice: ${slicePos.toFixed(1)} a₀`;

    if (bondLength > 0) {
        const systemNames = {
            h2_ion: 'H₂⁺',
            h2_molecule: 'H₂'
        };
        info += ` | ${systemNames[system]} bond length: ${bondLength.toFixed(2)} a₀`;
    }

    elements.densityInfo.textContent = info;
}

// Clear density canvas
function clearDensity() {
    const canvas = elements.densityCanvas;
    const ctx = canvas.getContext('2d');
    const colors = COLOR_MAPS[currentColorTheme] || COLOR_MAPS.viridis;
    const bgColor = colors[0];
    ctx.fillStyle = `rgb(${bgColor[0]}, ${bgColor[1]}, ${bgColor[2]})`;
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    elements.densityInfo.textContent = '';
}

// Draw nuclear positions and bond on density plot
function drawNucleiAndBond(ctx, plane, slicePos, bounds, size) {
    const system = elements.system.value;
    const bondLength = BOND_LENGTHS[system];
    let nuclei = [];

    if (system === 'hydrogen') {
        nuclei = [[0, 0, 0]];
    } else {
        const r = bondLength / 2;
        nuclei = [[-r, 0, 0], [r, 0, 0]];
    }

    // Project nuclei onto plane
    const projectedNuclei = [];
    for (const nuc of nuclei) {
        let u, v, w;
        if (plane === 'xy') {
            [u, v, w] = [nuc[0], nuc[1], nuc[2]];
        } else if (plane === 'xz') {
            [u, v, w] = [nuc[0], nuc[2], nuc[1]];
        } else {
            [u, v, w] = [nuc[1], nuc[2], nuc[0]];
        }

        // Check if nucleus is near the slice
        const nearSlice = Math.abs(w - slicePos) <= 0.5;

        const x = ((u + bounds) / (2 * bounds)) * size;
        const y = ((v + bounds) / (2 * bounds)) * size;

        projectedNuclei.push({ x, y, nearSlice, u, v });
    }

    // Draw bond line if both nuclei are visible and it's a 2-atom system
    if (nuclei.length === 2 && projectedNuclei[0].nearSlice && projectedNuclei[1].nearSlice) {
        ctx.strokeStyle = 'rgba(255, 255, 255, 0.7)';
        ctx.lineWidth = 2;
        ctx.setLineDash([5, 3]);
        ctx.beginPath();
        ctx.moveTo(projectedNuclei[0].x, projectedNuclei[0].y);
        ctx.lineTo(projectedNuclei[1].x, projectedNuclei[1].y);
        ctx.stroke();
        ctx.setLineDash([]);

        // Draw bond length label
        const midX = (projectedNuclei[0].x + projectedNuclei[1].x) / 2;
        const midY = (projectedNuclei[0].y + projectedNuclei[1].y) / 2;

        ctx.fillStyle = 'rgba(255, 255, 255, 0.9)';
        ctx.font = '12px monospace';
        ctx.textAlign = 'center';
        ctx.fillText(`R = ${bondLength.toFixed(2)} a₀`, midX, midY - 10);
    }

    // Draw nuclei
    for (const nuc of projectedNuclei) {
        if (!nuc.nearSlice) continue;

        ctx.fillStyle = '#fff';
        ctx.strokeStyle = '#fff';
        ctx.lineWidth = 2;

        ctx.beginPath();
        ctx.arc(nuc.x, nuc.y, 6, 0, 2 * Math.PI);
        ctx.fill();
        ctx.stroke();

        // Draw + symbol for nucleus
        ctx.fillStyle = '#000';
        ctx.font = 'bold 10px sans-serif';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillText('+', nuc.x, nuc.y);
    }
}

// Start the application
main().catch(console.error);
