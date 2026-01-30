//! Markdown to HTML conversion with KaTeX math rendering.

use pulldown_cmark::{CodeBlockKind, Event, Options, Parser as MdParser, Tag, TagEnd};
use std::collections::HashMap;

/// Tracks section numbers for hierarchical heading numbering.
struct SectionCounter {
    h2: usize,
    h3: usize,
    h4: usize,
}

impl SectionCounter {
    fn new() -> Self {
        Self { h2: 0, h3: 0, h4: 0 }
    }

    /// Increment counter for the given heading level and return the section number string.
    fn increment(&mut self, level: u8) -> String {
        match level {
            2 => {
                self.h2 += 1;
                self.h3 = 0;
                self.h4 = 0;
                format!("{}.", self.h2)
            }
            3 => {
                self.h3 += 1;
                self.h4 = 0;
                format!("{}.{}.", self.h2, self.h3)
            }
            4 => {
                self.h4 += 1;
                format!("{}.{}.{}.", self.h2, self.h3, self.h4)
            }
            _ => String::new(),
        }
    }
}

/// Convert markdown content to HTML with math rendering.
pub fn markdown_to_html(content: &str) -> Result<String, Box<dyn std::error::Error>> {
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
    let mut section_counter = SectionCounter::new(); // Track section numbers
    let mut slug_to_number: HashMap<String, String> = HashMap::new(); // Map slug to section number for TOC

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
                                // Also add original LaTeX to heading_text for slug generation
                                if let Some(original_latex) = inline_math_map.get(&id) {
                                    heading_text.push_str(original_latex);
                                }
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
            Event::Start(Tag::Image { .. }) => {
                // Pass through image tags as-is
                output_events.push(event.clone());
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

                // Generate section number (skip h1 and "Table of Contents")
                let section_num = if heading_level >= 2
                    && heading_level <= 4
                    && heading_text.trim() != "Table of Contents"
                {
                    section_counter.increment(heading_level)
                } else {
                    String::new()
                };

                // Store mapping for TOC post-processing
                if !section_num.is_empty() {
                    slug_to_number.insert(slug.clone(), section_num.clone());
                }

                let heading_html = if section_num.is_empty() {
                    format!(
                        "<h{} id=\"{}\">{}</h{}>",
                        heading_level, slug, heading_content, heading_level
                    )
                } else {
                    format!(
                        "<h{} id=\"{}\"><a href=\"#table-of-contents\" class=\"section-num\">{}</a> {}</h{}>",
                        heading_level, slug, section_num, heading_content, heading_level
                    )
                };
                output_events.push(Event::Html(heading_html.into()));
                continue;
            }
            _ => output_events.push(event.clone()),
        }
    }

    let mut html_output = String::new();
    pulldown_cmark::html::push_html(&mut html_output, output_events.into_iter());

    // Post-process to add section numbers to TOC links
    let html_output = add_numbers_to_toc(&html_output, &slug_to_number);

    Ok(html_output)
}

/// Add section numbers to TOC links.
fn add_numbers_to_toc(html: &str, slug_to_number: &HashMap<String, String>) -> String {
    let mut result = html.to_string();
    for (slug, number) in slug_to_number {
        // Try both unencoded and URL-encoded versions of the slug
        let patterns = [
            format!("<a href=\"#{}\">", slug),
            format!("<a href=\"#{}\">", url_encode_slug(slug)),
        ];

        for pattern in &patterns {
            if let Some(pos) = result.find(pattern) {
                let insert_pos = pos + pattern.len();
                let number_span = format!("<span class=\"toc-num\">{}</span> ", number);
                result.insert_str(insert_pos, &number_span);
                break; // Only insert once per slug
            }
        }
    }
    result
}

/// URL-encode special characters in slugs (for matching TOC links).
fn url_encode_slug(slug: &str) -> String {
    slug.chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '-' || c == '_' {
                c.to_string()
            } else {
                // Encode non-ASCII characters
                let mut buf = [0u8; 4];
                let encoded = c.encode_utf8(&mut buf);
                encoded.bytes().map(|b| format!("%{:02X}", b)).collect()
            }
        })
        .collect()
}

/// Extract inline math from GitHub syntax $`...`$ and replace with placeholders.
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

/// Render LaTeX math using KaTeX.
fn render_katex(latex: &str, display_mode: bool) -> Result<String, Box<dyn std::error::Error>> {
    let opts = katex::Opts::builder().display_mode(display_mode).build()?;
    let rendered = katex::render_with_opts(latex, &opts)?;

    if display_mode {
        Ok(format!("<div class=\"math-display\">{}</div>", rendered))
    } else {
        Ok(format!("<span class=\"math-inline\">{}</span>", rendered))
    }
}

/// Simple HTML escaping.
pub fn html_escape(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

/// Generate a URL-safe slug from heading text (matches GitHub's behavior).
fn slugify(text: &str) -> String {
    text.to_lowercase()
        .chars()
        .map(|c| match c {
            // Keep alphanumeric, accented chars, and underscores
            'a'..='z' | '0'..='9' | '_' => c,
            'ö' | 'ò' | 'ó' | 'ô' | 'õ' => c,
            'ä' | 'à' | 'á' | 'â' | 'ã' => c,
            'ü' | 'ù' | 'ú' | 'û' => c,
            'ë' | 'è' | 'é' | 'ê' => c,
            'ï' | 'ì' | 'í' | 'î' => c,
            'ñ' | 'ç' => c,
            // Convert spaces and hyphens to hyphens
            ' ' | '-' => '-',
            // Remove apostrophes entirely
            '\'' | '\u{2018}' | '\u{2019}' => '\0',
            // Remove other characters
            _ => '\0',
        })
        .filter(|&c| c != '\0')
        .collect::<String>()
        .split('-')
        .filter(|s| !s.is_empty())
        .collect::<Vec<_>>()
        .join("-")
}
