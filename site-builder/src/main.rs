//! Site builder for Diffusion Monte Carlo documentation and web simulator.
//!
//! # Commands
//!
//! ```bash
//! # Build WASM module
//! cargo run -p site-builder -- wasm
//!
//! # Build the complete site
//! cargo run -p site-builder -- build
//!
//! # Build and serve with dev server
//! cargo run -p site-builder -- serve
//!
//! # Serve on custom port
//! cargo run -p site-builder -- serve --port 3000
//! ```

use clap::{Parser, Subcommand};
use maud::{html, Markup, PreEscaped, DOCTYPE};
use pulldown_cmark::{CodeBlockKind, Event, Options, Parser as MdParser, Tag, TagEnd};
use std::collections::HashMap;
use std::fs;
use std::io::Read;
use std::path::Path;
use std::process::Command;
use tiny_http::{Response, Server};

#[derive(Parser)]
#[command(name = "site-builder")]
#[command(about = "Build and serve the Diffusion MC website")]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    /// Build WASM module using wasm-pack
    Wasm,
    /// Build the complete site (documentation + simulator)
    Build,
    /// Build and serve the site with a development server
    Serve {
        /// Port to serve on
        #[arg(short, long, default_value = "8080")]
        port: u16,
    },
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    match cli.command {
        Some(Commands::Wasm) => {
            build_wasm()?;
        }
        Some(Commands::Build) => {
            build_site()?;
        }
        Some(Commands::Serve { port }) => {
            build_site()?;
            serve_site(port)?;
        }
        None => {
            // Default: just build the site
            build_site()?;
        }
    }

    Ok(())
}

/// Build WASM module using wasm-pack
fn build_wasm() -> Result<(), Box<dyn std::error::Error>> {
    println!("Building WASM module...");

    // Check if wasm-pack is installed
    let wasm_pack_check = Command::new("wasm-pack").arg("--version").output();

    if wasm_pack_check.is_err() {
        println!("wasm-pack not found. Installing...");
        let status = Command::new("cargo")
            .args(["install", "wasm-pack"])
            .status()?;

        if !status.success() {
            return Err("Failed to install wasm-pack".into());
        }
    }

    // Build the WASM module
    let status = Command::new("wasm-pack")
        .args(["build", "--target", "web", "--out-dir", "www/pkg"])
        .current_dir("dmc-web")
        .status()?;

    if !status.success() {
        return Err("wasm-pack build failed".into());
    }

    println!("WASM build complete! Output in dmc-web/www/pkg/");
    Ok(())
}

/// Build the complete site
fn build_site() -> Result<(), Box<dyn std::error::Error>> {
    // Read the README.md file
    let readme_path = Path::new("README.md");
    let readme_content = fs::read_to_string(readme_path)?;

    // Convert markdown to HTML with math rendering
    let body_html = markdown_to_html(&readme_content)?;

    // Generate the full page using maud
    let page = generate_page(&body_html);

    // Create dist directory and write output
    let dist_path = Path::new("dist");
    fs::create_dir_all(dist_path)?;
    fs::write(dist_path.join("index.html"), page.into_string())?;

    // Copy images directory if it exists
    let images_src = Path::new("images");
    let images_dest = dist_path.join("images");
    if images_src.exists() {
        fs::create_dir_all(&images_dest)?;
        for entry in fs::read_dir(images_src)? {
            let entry = entry?;
            let path = entry.path();
            if path.is_file() {
                let filename = path.file_name().unwrap();
                fs::copy(&path, images_dest.join(filename))?;
            }
        }
    }

    // Copy simulator files from dmc-web/www
    let simulator_src = Path::new("dmc-web/www");
    let simulator_dest = dist_path.join("simulator");

    // Check if WASM has been built
    let pkg_path = simulator_src.join("pkg");
    if !pkg_path.exists() {
        println!("WASM not built yet. Building now...");
        build_wasm()?;
    }

    if simulator_src.exists() {
        copy_simulator_files(simulator_src, &simulator_dest)?;
        println!("Simulator files copied to dist/simulator/");
    }

    println!("Site built successfully! Output in dist/");
    Ok(())
}

/// Serve the site with a simple HTTP server
fn serve_site(port: u16) -> Result<(), Box<dyn std::error::Error>> {
    let addr = format!("0.0.0.0:{}", port);
    let server = Server::http(&addr).map_err(|e| format!("Failed to start server: {}", e))?;

    println!("Serving site at http://localhost:{}", port);
    println!("Press Ctrl+C to stop");

    for request in server.incoming_requests() {
        let url = request.url().to_string();
        let path = if url == "/" {
            "dist/index.html".to_string()
        } else if url.ends_with('/') {
            format!("dist{}index.html", url)
        } else {
            format!("dist{}", url)
        };

        let response = match serve_file(&path) {
            Ok(resp) => resp,
            Err(_) => {
                // Try with .html extension
                match serve_file(&format!("{}.html", path)) {
                    Ok(resp) => resp,
                    Err(_) => {
                        let not_found = "404 Not Found";
                        Response::from_string(not_found).with_status_code(404)
                    }
                }
            }
        };

        let _ = request.respond(response);
    }

    Ok(())
}

/// Serve a single file
fn serve_file(path: &str) -> Result<Response<std::io::Cursor<Vec<u8>>>, std::io::Error> {
    let mut file = fs::File::open(path)?;
    let mut content = Vec::new();
    file.read_to_end(&mut content)?;

    let content_type = match Path::new(path).extension().and_then(|e| e.to_str()) {
        Some("html") => "text/html; charset=utf-8",
        Some("css") => "text/css; charset=utf-8",
        Some("js") => "application/javascript; charset=utf-8",
        Some("wasm") => "application/wasm",
        Some("json") => "application/json; charset=utf-8",
        Some("png") => "image/png",
        Some("jpg") | Some("jpeg") => "image/jpeg",
        Some("svg") => "image/svg+xml",
        _ => "application/octet-stream",
    };

    Ok(Response::from_data(content)
        .with_header(tiny_http::Header::from_bytes("Content-Type", content_type).unwrap())
        .with_header(
            tiny_http::Header::from_bytes("Cross-Origin-Opener-Policy", "same-origin").unwrap(),
        )
        .with_header(
            tiny_http::Header::from_bytes("Cross-Origin-Embedder-Policy", "require-corp").unwrap(),
        ))
}

/// Copy simulator files (HTML, CSS, JS, WASM) to dist
fn copy_simulator_files(src: &Path, dest: &Path) -> Result<(), Box<dyn std::error::Error>> {
    fs::create_dir_all(dest)?;

    // Copy main web files
    for file in ["index.html", "style.css", "main.js"] {
        let src_file = src.join(file);
        if src_file.exists() {
            fs::copy(&src_file, dest.join(file))?;
        }
    }

    // Copy pkg directory (WASM and generated JS)
    let pkg_src = src.join("pkg");
    let pkg_dest = dest.join("pkg");
    if pkg_src.exists() {
        fs::create_dir_all(&pkg_dest)?;
        copy_dir_contents(&pkg_src, &pkg_dest)?;
    }

    Ok(())
}

/// Recursively copy directory contents
fn copy_dir_contents(src: &Path, dest: &Path) -> Result<(), Box<dyn std::error::Error>> {
    for entry in fs::read_dir(src)? {
        let entry = entry?;
        let path = entry.path();
        let filename = path.file_name().unwrap();
        let dest_path = dest.join(filename);

        if path.is_dir() {
            fs::create_dir_all(&dest_path)?;
            copy_dir_contents(&path, &dest_path)?;
        } else {
            fs::copy(&path, &dest_path)?;
        }
    }
    Ok(())
}

/// Convert markdown to HTML, rendering math with KaTeX
fn markdown_to_html(content: &str) -> Result<String, Box<dyn std::error::Error>> {
    // Step 1: Extract all inline math and replace with numbered placeholders
    let (processed, inline_math_map) = extract_inline_math(content);

    // Step 2: Pre-render all inline math with KaTeX
    let mut rendered_math: HashMap<usize, String> = HashMap::new();
    for (id, latex) in &inline_math_map {
        match render_katex(latex, false) {
            Ok(html) => {
                rendered_math.insert(*id, html);
            }
            Err(e) => {
                eprintln!("KaTeX error for inline math '{}': {}", latex, e);
                rendered_math.insert(
                    *id,
                    format!("<code class=\"math-error\">{}</code>", html_escape(latex)),
                );
            }
        }
    }

    // Step 3: Parse markdown with pulldown-cmark
    let options = Options::ENABLE_TABLES | Options::ENABLE_STRIKETHROUGH | Options::ENABLE_TASKLISTS;

    let parser = MdParser::new_ext(&processed, options);

    // Process events, handling math code blocks and headings
    let mut in_math_block = false;
    let mut math_content = String::new();
    let mut in_heading = false;
    let mut heading_level: u8 = 0;
    let mut heading_content = String::new(); // Accumulates HTML content for heading
    let mut heading_text = String::new(); // Plain text for slug generation
    let mut slug_counts: HashMap<String, usize> = HashMap::new(); // Track duplicate headings
    let mut output_events: Vec<Event> = Vec::new();

    for event in parser {
        match &event {
            Event::Start(Tag::CodeBlock(CodeBlockKind::Fenced(lang))) => {
                if lang.as_ref() == "math" {
                    in_math_block = true;
                    math_content.clear();
                    continue;
                }
                output_events.push(event);
            }
            Event::End(TagEnd::CodeBlock) => {
                if in_math_block {
                    in_math_block = false;
                    match render_katex(&math_content, true) {
                        Ok(rendered) => {
                            output_events.push(Event::Html(rendered.into()));
                        }
                        Err(e) => {
                            eprintln!("KaTeX error for display math: {}", e);
                            let fallback =
                                format!("<pre class=\"math-error\">{}</pre>", html_escape(&math_content));
                            output_events.push(Event::Html(fallback.into()));
                        }
                    }
                    continue;
                }
                output_events.push(event);
            }
            Event::Text(text) => {
                if in_math_block {
                    math_content.push_str(text);
                    continue;
                }
                if in_heading {
                    heading_text.push_str(text);
                    heading_content.push_str(&html_escape(text));
                    continue;
                }
                output_events.push(event.clone());
            }
            Event::Code(code) => {
                if code.starts_with("INLINEMATH") {
                    if let Ok(id) = code.trim_start_matches("INLINEMATH").parse::<usize>() {
                        if let Some(rendered) = rendered_math.get(&id) {
                            if in_heading {
                                // Add rendered math inline within the heading
                                heading_content.push_str(rendered);
                            } else {
                                output_events.push(Event::Html(rendered.clone().into()));
                            }
                            continue;
                        }
                    }
                }
                if in_heading {
                    // Regular inline code in heading
                    heading_content.push_str("<code>");
                    heading_content.push_str(&html_escape(code));
                    heading_content.push_str("</code>");
                    heading_text.push_str(code);
                    continue;
                }
                output_events.push(event.clone());
            }
            Event::Start(Tag::Image {
                link_type,
                dest_url,
                title,
                id,
            }) => {
                let new_dest = if dest_url.contains("raw.githubusercontent.com") {
                    if let Some(filename) = dest_url.split('/').last() {
                        format!("images/{}", filename).into()
                    } else {
                        dest_url.clone()
                    }
                } else {
                    dest_url.clone()
                };
                output_events.push(Event::Start(Tag::Image {
                    link_type: *link_type,
                    dest_url: new_dest,
                    title: title.clone(),
                    id: id.clone(),
                }));
            }
            Event::Start(Tag::Heading { level, .. }) => {
                in_heading = true;
                heading_level = *level as u8;
                heading_text.clear();
                heading_content.clear();
                continue;
            }
            Event::End(TagEnd::Heading(_)) => {
                in_heading = false;
                let base_slug = slugify(&heading_text);
                // Handle duplicate headings by appending -1, -2, etc.
                let count = slug_counts.entry(base_slug.clone()).or_insert(0);
                let slug = if *count == 0 {
                    base_slug
                } else {
                    format!("{}-{}", base_slug, count)
                };
                *slug_counts.get_mut(&slugify(&heading_text)).unwrap() += 1;
                let heading_html = format!(
                    "<h{} id=\"{}\">{}</h{}>",
                    heading_level,
                    slug,
                    heading_content,
                    heading_level
                );
                output_events.push(Event::Html(heading_html.into()));
                continue;
            }
            _ => output_events.push(event.clone()),
        }
    }

    let mut html_output = String::new();
    pulldown_cmark::html::push_html(&mut html_output, output_events.into_iter());

    Ok(html_output)
}

/// Extract inline math from GitHub syntax $`...`$ and replace with placeholders
fn extract_inline_math(content: &str) -> (String, HashMap<usize, String>) {
    let mut result = String::with_capacity(content.len());
    let mut math_map: HashMap<usize, String> = HashMap::new();
    let mut math_id = 0;
    let mut chars = content.chars().peekable();

    while let Some(c) = chars.next() {
        if c == '$' && chars.peek() == Some(&'`') {
            chars.next();
            let mut math = String::new();

            loop {
                match chars.next() {
                    Some('`') if chars.peek() == Some(&'$') => {
                        chars.next();
                        break;
                    }
                    Some(mc) => math.push(mc),
                    None => break,
                }
            }

            math_map.insert(math_id, math);
            result.push('`');
            result.push_str(&format!("INLINEMATH{}", math_id));
            result.push('`');
            math_id += 1;
        } else {
            result.push(c);
        }
    }

    (result, math_map)
}

/// Render LaTeX math using KaTeX
fn render_katex(latex: &str, display_mode: bool) -> Result<String, Box<dyn std::error::Error>> {
    let opts = katex::Opts::builder().display_mode(display_mode).build()?;
    let rendered = katex::render_with_opts(latex, &opts)?;

    if display_mode {
        Ok(format!("<div class=\"math-display\">{}</div>", rendered))
    } else {
        Ok(format!("<span class=\"math-inline\">{}</span>", rendered))
    }
}

/// Simple HTML escaping
fn html_escape(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

/// Generate a URL-safe slug from heading text
fn slugify(text: &str) -> String {
    text.to_lowercase()
        .chars()
        .map(|c| match c {
            'a'..='z' | '0'..='9' => c,
            ' ' | '-' | '_' => '-',
            // Convert common accented characters to ASCII equivalents
            'ö' | 'ò' | 'ó' | 'ô' | 'õ' => 'o',
            'ä' | 'à' | 'á' | 'â' | 'ã' => 'a',
            'ü' | 'ù' | 'ú' | 'û' => 'u',
            'ë' | 'è' | 'é' | 'ê' => 'e',
            'ï' | 'ì' | 'í' | 'î' => 'i',
            'ñ' => 'n',
            'ç' => 'c',
            // Remove apostrophes entirely (no hyphen)
            '\'' | '\u{2018}' | '\u{2019}' => '\0',
            _ => '-',
        })
        .filter(|&c| c != '\0')
        .collect::<String>()
        .split('-')
        .filter(|s| !s.is_empty())
        .collect::<Vec<_>>()
        .join("-")
}

/// Generate the complete HTML page
fn generate_page(body_content: &str) -> Markup {
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
