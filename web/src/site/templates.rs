//! Maud templates for HTML page generation.

use maud::{html, Markup, PreEscaped, DOCTYPE};
use super::assets::SIMULATOR_CSS;
use super::markdown::markdown_to_html;

/// Generate the documentation page from README markdown.
pub fn generate_documentation_page(readme_content: &str) -> Result<String, Box<dyn std::error::Error>> {
    let body_html = markdown_to_html(readme_content)?;
    Ok(documentation_page(&body_html).into_string())
}

/// Generate the simulator page.
pub fn generate_simulator_page() -> String {
    simulator_page().into_string()
}

/// Documentation page template.
fn documentation_page(body_content: &str) -> Markup {
    html! {
        (DOCTYPE)
        html lang="en" data-theme="light" {
            head {
                meta charset="UTF-8";
                meta name="viewport" content="width=device-width, initial-scale=1.0";
                title { "Diffusion Monte Carlo" }
                link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/@picocss/pico@2/css/pico.classless.min.css";
                link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.7/dist/katex.min.css";
                style {
                    ":root { --pico-font-size: 106.25%; }"
                    ".math-display { overflow-x: auto; padding: 1rem 0; text-align: center; }"
                    ".math-inline { white-space: nowrap; }"
                    ".math-error { color: var(--pico-del-color); }"
                    "main { padding-top: 2rem; }"
                    "nav { position: sticky; top: 0; background: var(--pico-background-color); z-index: 100; border-bottom: 1px solid var(--pico-muted-border-color); padding: 0.75rem 0; }"
                    "nav .container { width: 100%; max-width: 1200px; margin: 0 auto; padding: 0 2rem; display: flex; justify-content: space-between; align-items: center; box-sizing: border-box; }"
                    "nav .brand { font-weight: 600; font-size: 1.1rem; color: var(--pico-primary); text-decoration: none; white-space: nowrap; }"
                    "nav .brand:hover { color: var(--pico-primary-hover); }"
                    "nav ul { display: flex; gap: 2rem; margin: 0; padding: 0; list-style: none; flex-shrink: 0; }"
                    "nav li { margin: 0; }"
                    "nav a { text-decoration: none; color: var(--pico-muted-color); padding: 0.5rem 0.75rem; transition: color 0.2s; }"
                    "nav a:hover { color: var(--pico-color); }"
                    "nav a.active { color: var(--pico-color); font-weight: 600; }"
                }
            }
            body {
                nav {
                    div class="container" {
                        a class="brand" href="index.html" { "Diffusion MC" }
                        ul {
                            li { a class="active" href="index.html" { "Documentation" } }
                            li { a href="simulator/index.html" { "Simulator" } }
                        }
                    }
                }
                main {
                    (PreEscaped(body_content))
                }
                footer {
                    p { "Copyright © 2017 Subeen Pang" }
                }
            }
        }
    }
}

/// Simulator page template.
fn simulator_page() -> Markup {
    html! {
        (DOCTYPE)
        html lang="en" {
            head {
                meta charset="UTF-8";
                meta name="viewport" content="width=device-width, initial-scale=1.0";
                title { "Diffusion Monte Carlo Simulator" }
                link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/@picocss/pico@2/css/pico.min.css";
                script src="https://cdn.jsdelivr.net/npm/chart.js@4" {}
                style { (SIMULATOR_CSS) }
            }
            body {
                (nav_component("../index.html", false))
                main class="container" {
                    header {
                        h1 { "Diffusion Monte Carlo Simulator" }
                        p { "Real-time quantum ground state calculations in the browser" }
                    }

                    div class="grid" {
                        (control_panel())
                        (statistics_panel())
                    }

                    div class="grid" {
                        section {
                            article {
                                header { h3 { "Energy vs. Step" } }
                                canvas id="energy_chart" {}
                            }
                        }
                        section {
                            article {
                                header { h3 { "Population vs. Step" } }
                                canvas id="population_chart" {}
                            }
                        }
                    }

                    (density_section())

                    footer {
                        p {
                            small {
                                "Powered by "
                                a href="https://github.com/subinbg/Diffusion-MC" { "Diffusion-MC" }
                            }
                        }
                    }
                }
                (noscript_fallback())
                (module_error_handler())
                script type="module" src="main.js" onerror="showServerRequiredMessage()" {}
            }
        }
    }
}

/// Navigation component.
fn nav_component(docs_href: &str, docs_active: bool) -> Markup {
    html! {
        nav {
            div class="nav-container" {
                a class="brand" href=(docs_href) { "Diffusion MC" }
                ul {
                    li {
                        a class=[docs_active.then_some("active")] href=(docs_href) { "Documentation" }
                    }
                    li {
                        a class=[(!docs_active).then_some("active")]
                          href=[if docs_active { Some("simulator/index.html") } else { Some("index.html") }]
                        { "Simulator" }
                    }
                }
            }
        }
    }
}

/// Control panel with configuration inputs.
fn control_panel() -> Markup {
    html! {
        section id="controls" {
            article {
                header { h3 { "Configuration" } }

                label for="system" { "Quantum System" }
                select id="system" name="system" {
                    option value="hydrogen" { "Hydrogen (H)" }
                    option value="h2_ion" { "H2+ Ion" }
                    option value="h2_molecule" { "H2 Molecule" }
                }

                label for="algorithm" { "Algorithm" }
                select id="algorithm" name="algorithm" {
                    option value="importance_sampled" selected { "Importance Sampled" }
                    option value="pure" { "Pure DMC" }
                }

                div id="zero_variance_warning" class="warning-box" style="display: none;" {
                    strong { "Note: " }
                    "With the exact trial wavefunction (α = 1), the hydrogen ground state has zero variance."
                }

                div class="grid" {
                    div {
                        label for="num_walkers" { "Walkers" }
                        input type="number" id="num_walkers" value="1000" min="100" max="10000" step="100";
                    }
                    div {
                        label for="time_step" { "Time Step" }
                        input type="number" id="time_step" value="0.01" min="0.001" max="0.1" step="0.001";
                    }
                }

                div class="grid" {
                    div {
                        label for="total_steps" { "Total Steps" }
                        input type="number" id="total_steps" value="5000" min="1000" max="100000" step="1000";
                    }
                    div {
                        label for="equilibration_steps" { "Equilibration" }
                        input type="number" id="equilibration_steps" value="500" min="0" max="10000" step="100";
                    }
                }

                label for="update_interval" {
                    "Update Interval: "
                    span id="update_interval_value" { "10" }
                    " steps"
                }
                input type="range" id="update_interval" min="1" max="100" value="10";

                footer {
                    div class="grid" {
                        button id="start_btn" class="primary" { "Start" }
                        button id="pause_btn" class="secondary" disabled { "Pause" }
                        button id="reset_btn" class="contrast" { "Reset" }
                    }
                }
            }
        }
    }
}

/// Statistics panel with results display.
fn statistics_panel() -> Markup {
    html! {
        section id="statistics" {
            article {
                header { h3 { "Results" } }

                div class="stats-grid" {
                    div class="stat-item" {
                        span class="stat-label" { "Step" }
                        span class="stat-value" id="stat_step" { "0" }
                        span class="stat-unit" { "/ " span id="stat_total" { "5000" } }
                    }
                    div class="stat-item" {
                        span class="stat-label" { "Phase" }
                        span class="stat-value" id="stat_phase" { "Ready" }
                    }
                    div class="stat-item" {
                        span class="stat-label" { "Energy" }
                        span class="stat-value" id="stat_energy" { "-" }
                        span class="stat-unit" { "Ha" }
                    }
                    div class="stat-item" {
                        span class="stat-label" { "Error" }
                        span class="stat-value" id="stat_error" { "-" }
                        span class="stat-unit" { "Ha" }
                    }
                    div class="stat-item" {
                        span class="stat-label" { "Reference" }
                        span class="stat-value" id="stat_exact" { "-0.5" }
                        span class="stat-unit" { "Ha" }
                    }
                    div class="stat-item" {
                        span class="stat-label" { "Population" }
                        span class="stat-value" id="stat_population" { "0" }
                        span class="stat-unit" { "walkers" }
                    }
                    div class="stat-item" {
                        span class="stat-label" { "Acceptance" }
                        span class="stat-value" id="stat_acceptance" { "-" }
                        span class="stat-unit" { "%" }
                    }
                    div class="stat-item" {
                        span class="stat-label" { "Deviation" }
                        span class="stat-value" id="stat_deviation" { "-" }
                        span class="stat-unit" { "%" }
                    }
                }

                progress id="progress" value="0" max="100" {}
            }
        }
    }
}

/// Density visualization section.
fn density_section() -> Markup {
    html! {
        section {
            article {
                header {
                    h3 { "Electron Density |ψ|²" }
                    div class="density-controls" {
                        label for="density_plane" { "Plane:" }
                        select id="density_plane" {
                            option value="xy" selected { "XY" }
                            option value="xz" { "XZ" }
                            option value="yz" { "YZ" }
                        }
                        label for="slice_position" {
                            "Slice: "
                            span id="slice_value" { "0.0" }
                            " a₀"
                        }
                        input type="range" id="slice_position" min="-3" max="3" value="0" step="0.1";
                        label for="color_theme" { "Color:" }
                        select id="color_theme" {
                            option value="magma" selected { "Magma" }
                            option value="viridis" { "Viridis" }
                            option value="plasma" { "Plasma" }
                            option value="inferno" { "Inferno" }
                            option value="grayscale" { "Grayscale" }
                        }
                    }
                }
                div class="density-container" {
                    canvas id="density_canvas" width="400" height="400" {}
                    div class="colorbar-container" {
                        span class="colorbar-label" id="colorbar_max" { "1.0" }
                        div id="colorbar" {}
                        span class="colorbar-label" id="colorbar_min" { "0.0" }
                    }
                }
                div class="density-info" id="density_info" {}
            }
        }
    }
}

/// Noscript fallback message.
fn noscript_fallback() -> Markup {
    html! {
        noscript {
            article style="margin: 2rem;" {
                h2 { "JavaScript Required" }
                p { "This simulator requires JavaScript to run. Please enable JavaScript in your browser settings." }
            }
        }
    }
}

/// Module error handler script.
fn module_error_handler() -> Markup {
    html! {
        script {
            (PreEscaped(r#"
                window.addEventListener('error', function(e) {
                    if (e.message && e.message.includes('module')) {
                        showServerRequiredMessage();
                    }
                });

                function showServerRequiredMessage() {
                    document.body.innerHTML = `
                        <main class="container" style="margin-top: 2rem;">
                            <article>
                                <header><h2>Server Required</h2></header>
                                <p>This simulator uses WebAssembly and ES modules, which require a web server to run properly.</p>
                                <p>Opening the HTML file directly from your filesystem won't work due to browser security restrictions.</p>
                                <h3>To run the simulator:</h3>
                                <pre><code>cd Diffusion-MC
cargo run -p web --features cli -- serve</code></pre>
                                <p>Then open <a href="http://localhost:8080/simulator/">http://localhost:8080/simulator/</a></p>
                            </article>
                        </main>
                    `;
                }
            "#))
        }
    }
}
