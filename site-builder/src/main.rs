use maud::{html, Markup, PreEscaped, DOCTYPE};
use pulldown_cmark::{CodeBlockKind, Event, Options, Parser, Tag, TagEnd};
use std::collections::HashMap;
use std::fs;
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
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

    println!("Site built successfully! Output in dist/");
    Ok(())
}

/// Convert markdown to HTML, rendering math with KaTeX
fn markdown_to_html(content: &str) -> Result<String, Box<dyn std::error::Error>> {
    // Step 1: Extract all inline math and replace with numbered placeholders
    // This prevents markdown parser from interfering with math content
    let (processed, inline_math_map) = extract_inline_math(content);

    // Step 2: Pre-render all inline math with KaTeX
    let mut rendered_math: HashMap<usize, String> = HashMap::new();
    for (id, latex) in &inline_math_map {
        match render_katex(latex, false) {
            Ok(html) => { rendered_math.insert(*id, html); }
            Err(e) => {
                eprintln!("KaTeX error for inline math '{}': {}", latex, e);
                // Fallback: show raw LaTeX in a code element
                rendered_math.insert(*id, format!("<code class=\"math-error\">{}</code>", html_escape(latex)));
            }
        }
    }

    // Step 3: Parse markdown with pulldown-cmark
    let options = Options::ENABLE_TABLES
        | Options::ENABLE_STRIKETHROUGH
        | Options::ENABLE_TASKLISTS;

    let parser = Parser::new_ext(&processed, options);

    // Process events, handling math code blocks
    let mut in_math_block = false;
    let mut math_content = String::new();
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
                    // Render display math block with KaTeX
                    match render_katex(&math_content, true) {
                        Ok(rendered) => {
                            output_events.push(Event::Html(rendered.into()));
                        }
                        Err(e) => {
                            eprintln!("KaTeX error for display math: {}", e);
                            let fallback = format!("<pre class=\"math-error\">{}</pre>", html_escape(&math_content));
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
                output_events.push(event.clone());
            }
            Event::Code(code) => {
                // Check if this is one of our math placeholders
                if code.starts_with("INLINEMATH") {
                    if let Ok(id) = code.trim_start_matches("INLINEMATH").parse::<usize>() {
                        if let Some(rendered) = rendered_math.get(&id) {
                            output_events.push(Event::Html(rendered.clone().into()));
                            continue;
                        }
                    }
                }
                output_events.push(event.clone());
            }
            Event::Start(Tag::Image { link_type, dest_url, title, id }) => {
                // Fix image paths - convert raw GitHub URLs to local paths
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
            _ => output_events.push(event.clone()),
        }
    }

    // Convert events to HTML
    let mut html_output = String::new();
    pulldown_cmark::html::push_html(&mut html_output, output_events.into_iter());

    Ok(html_output)
}

/// Extract inline math from GitHub syntax $`...`$ and replace with placeholders
/// Returns the processed content and a map of placeholder IDs to math content
fn extract_inline_math(content: &str) -> (String, HashMap<usize, String>) {
    let mut result = String::with_capacity(content.len());
    let mut math_map: HashMap<usize, String> = HashMap::new();
    let mut math_id = 0;
    let mut chars = content.chars().peekable();

    while let Some(c) = chars.next() {
        if c == '$' && chars.peek() == Some(&'`') {
            chars.next(); // consume the `
            let mut math = String::new();

            // Collect until we find `$
            loop {
                match chars.next() {
                    Some('`') if chars.peek() == Some(&'$') => {
                        chars.next(); // consume the $
                        break;
                    }
                    Some(mc) => math.push(mc),
                    None => break, // EOF, incomplete math
                }
            }

            // Store the math content and insert a placeholder
            // Use backticks to make it inline code, which markdown won't process further
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
    let opts = katex::Opts::builder()
        .display_mode(display_mode)
        .build()?;

    let rendered = katex::render_with_opts(latex, &opts)?;

    if display_mode {
        Ok(format!("<div class=\"math-display\">{}</div>", rendered))
    } else {
        Ok(format!("<span class=\"math-inline\">{}</span>", rendered))
    }
}

/// Simple HTML escaping for text content
fn html_escape(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

/// Generate the complete HTML page using maud
fn generate_page(body_content: &str) -> Markup {
    html! {
        (DOCTYPE)
        html lang="en" {
            head {
                meta charset="UTF-8";
                meta name="viewport" content="width=device-width, initial-scale=1.0";
                title { "Diffusion Monte Carlo" }

                // Simple.css - classless CSS framework
                link rel="stylesheet" href="https://cdn.simplecss.org/simple.min.css";

                // KaTeX CSS for math styling
                link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.7/dist/katex.min.css";

                // Minimal custom styles for math display
                style {
                    ".math-display { overflow-x: auto; padding: 1rem 0; text-align: center; }"
                    ".math-inline { white-space: nowrap; }"
                    ".math-error { color: red; }"
                }
            }
            body {
                main {
                    (PreEscaped(body_content))
                }
                footer {
                    p { "Generated with a Rust-based static site builder using maud and KaTeX" }
                }
            }
        }
    }
}
