use pulldown_cmark::{html, Event, Options, Parser, Tag, TagEnd, CodeBlockKind};
use std::fs;
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Read the README.md file
    let readme_path = Path::new("README.md");
    let readme_content = fs::read_to_string(readme_path)?;

    // Pre-process: Convert GitHub math syntax to MathJax-compatible syntax
    let processed_content = preprocess_math(&readme_content);

    // Set up pulldown-cmark options for GitHub Flavored Markdown
    let options = Options::ENABLE_TABLES
        | Options::ENABLE_STRIKETHROUGH
        | Options::ENABLE_TASKLISTS;

    // Parse markdown with custom handling for math code blocks
    let parser = Parser::new_ext(&processed_content, options);

    // Custom event processing to handle math blocks and fix image paths
    let mut in_math_block = false;
    let mut math_content = String::new();

    let events: Vec<Event> = parser
        .filter_map(|event| {
            match &event {
                Event::Start(Tag::CodeBlock(CodeBlockKind::Fenced(lang))) => {
                    if lang.as_ref() == "math" {
                        in_math_block = true;
                        math_content.clear();
                        return None;
                    }
                    Some(event)
                }
                Event::End(TagEnd::CodeBlock) => {
                    if in_math_block {
                        in_math_block = false;
                        // Return the math block as raw HTML with MathJax delimiters
                        let html = format!(
                            "<div class=\"math-display\">\\[\n{}\n\\]</div>",
                            math_content.trim()
                        );
                        return Some(Event::Html(html.into()));
                    }
                    Some(event)
                }
                Event::Text(text) => {
                    if in_math_block {
                        math_content.push_str(text);
                        return None;
                    }
                    Some(event.clone())
                }
                Event::Start(Tag::Image { link_type, dest_url, title, id }) => {
                    // Fix image paths - convert raw GitHub URLs to local paths
                    let new_dest = if dest_url.contains("raw.githubusercontent.com") {
                        // Extract just the filename from GitHub raw URLs
                        if let Some(filename) = dest_url.split('/').last() {
                            format!("images/{}", filename).into()
                        } else {
                            dest_url.clone()
                        }
                    } else {
                        dest_url.clone()
                    };
                    Some(Event::Start(Tag::Image {
                        link_type: *link_type,
                        dest_url: new_dest,
                        title: title.clone(),
                        id: id.clone(),
                    }))
                }
                _ => Some(event.clone()),
            }
        })
        .collect();

    // Convert to HTML
    let mut html_output = String::new();
    html::push_html(&mut html_output, events.into_iter());

    // Post-process to replace math placeholders with MathJax delimiters
    let html_output = postprocess_math(&html_output);

    // Wrap in HTML template
    let full_html = generate_html_page(&html_output);

    // Create dist directory and write output
    let dist_path = Path::new("dist");
    fs::create_dir_all(dist_path)?;
    fs::write(dist_path.join("index.html"), full_html)?;

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

/// Pre-process GitHub-flavored math syntax to MathJax-compatible syntax
/// Uses unique placeholders that won't be interpreted by markdown
const MATH_INLINE_START: &str = "MATHINLINESTART";
const MATH_INLINE_END: &str = "MATHINLINEEND";

fn preprocess_math(content: &str) -> String {
    let mut result = String::with_capacity(content.len());
    let mut chars = content.chars().peekable();

    while let Some(c) = chars.next() {
        if c == '$' && chars.peek() == Some(&'`') {
            // Found $` - start of inline math
            chars.next(); // consume the `
            let mut math = String::new();

            // Collect until we find `$
            while let Some(mc) = chars.next() {
                if mc == '`' && chars.peek() == Some(&'$') {
                    chars.next(); // consume the $
                    break;
                }
                math.push(mc);
            }

            // Output with placeholder markers
            result.push_str(MATH_INLINE_START);
            result.push_str(&math);
            result.push_str(MATH_INLINE_END);
        } else {
            result.push(c);
        }
    }

    result
}

/// Post-process HTML to replace math placeholders with MathJax delimiters
fn postprocess_math(html: &str) -> String {
    html.replace(MATH_INLINE_START, "\\(")
        .replace(MATH_INLINE_END, "\\)")
}

/// Generate the complete HTML page with styling
fn generate_html_page(content: &str) -> String {
    format!(
        r#"<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Diffusion Monte Carlo</title>

    <!-- MathJax for LaTeX rendering -->
    <script>
        MathJax = {{
            tex: {{
                inlineMath: [['\\(', '\\)']],
                displayMath: [['\\[', '\\]']],
                processEscapes: true
            }},
            options: {{
                skipHtmlTags: ['script', 'noscript', 'style', 'textarea', 'pre', 'code']
            }}
        }};
    </script>
    <script src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js" async></script>

    <!-- Highlight.js for code syntax highlighting -->
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/styles/github.min.css">
    <script src="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/highlight.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/languages/fortran.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/languages/python.min.js"></script>
    <script>hljs.highlightAll();</script>

    <style>
        :root {{
            --primary-color: #2c3e50;
            --accent-color: #3498db;
            --text-color: #333;
            --bg-color: #fafafa;
            --code-bg: #f4f4f4;
            --border-color: #e1e4e8;
            --max-width: 900px;
        }}

        * {{
            box-sizing: border-box;
        }}

        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            line-height: 1.7;
            color: var(--text-color);
            background-color: var(--bg-color);
            margin: 0;
            padding: 0;
        }}

        .container {{
            max-width: var(--max-width);
            margin: 0 auto;
            padding: 2rem;
            background: white;
            min-height: 100vh;
            box-shadow: 0 0 20px rgba(0,0,0,0.05);
        }}

        h1 {{
            font-size: 2.5rem;
            color: var(--primary-color);
            border-bottom: 3px solid var(--accent-color);
            padding-bottom: 0.5rem;
            margin-top: 0;
        }}

        h2 {{
            font-size: 1.8rem;
            color: var(--primary-color);
            margin-top: 2.5rem;
            padding-bottom: 0.3rem;
            border-bottom: 1px solid var(--border-color);
        }}

        h3 {{
            font-size: 1.4rem;
            color: var(--primary-color);
            margin-top: 2rem;
        }}

        h4 {{
            font-size: 1.2rem;
            color: var(--primary-color);
            margin-top: 1.5rem;
        }}

        h5 {{
            font-size: 1.1rem;
            color: var(--primary-color);
            margin-top: 1.2rem;
        }}

        p {{
            margin: 1rem 0;
        }}

        a {{
            color: var(--accent-color);
            text-decoration: none;
        }}

        a:hover {{
            text-decoration: underline;
        }}

        img {{
            max-width: 100%;
            height: auto;
            display: block;
            margin: 1.5rem auto;
            border-radius: 4px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}

        pre {{
            background: var(--code-bg);
            padding: 1rem;
            border-radius: 6px;
            overflow-x: auto;
            border: 1px solid var(--border-color);
        }}

        code {{
            font-family: 'SFMono-Regular', Consolas, 'Liberation Mono', Menlo, monospace;
            font-size: 0.9em;
        }}

        :not(pre) > code {{
            background: var(--code-bg);
            padding: 0.2em 0.4em;
            border-radius: 3px;
        }}

        blockquote {{
            margin: 1.5rem 0;
            padding: 0.5rem 1rem;
            border-left: 4px solid var(--accent-color);
            background: #f8f9fa;
            color: #555;
        }}

        blockquote p {{
            margin: 0.5rem 0;
        }}

        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 1.5rem 0;
            font-size: 0.95rem;
        }}

        th, td {{
            padding: 0.75rem 1rem;
            text-align: left;
            border: 1px solid var(--border-color);
        }}

        th {{
            background: var(--primary-color);
            color: white;
            font-weight: 600;
        }}

        tr:nth-child(even) {{
            background: #f8f9fa;
        }}

        tr:hover {{
            background: #f0f4f8;
        }}

        hr {{
            border: none;
            border-top: 1px solid var(--border-color);
            margin: 2rem 0;
        }}

        ul, ol {{
            padding-left: 1.5rem;
        }}

        li {{
            margin: 0.5rem 0;
        }}

        .math-display {{
            overflow-x: auto;
            padding: 1rem 0;
            text-align: center;
        }}

        /* Navigation styling for TOC */
        ul ul {{
            margin-top: 0.3rem;
        }}

        /* Strong text in figure captions */
        strong {{
            font-weight: 600;
        }}

        /* Footer */
        .footer {{
            margin-top: 4rem;
            padding-top: 2rem;
            border-top: 1px solid var(--border-color);
            text-align: center;
            color: #666;
            font-size: 0.9rem;
        }}

        /* Responsive */
        @media (max-width: 768px) {{
            .container {{
                padding: 1rem;
            }}

            h1 {{
                font-size: 2rem;
            }}

            h2 {{
                font-size: 1.5rem;
            }}

            table {{
                font-size: 0.85rem;
            }}

            th, td {{
                padding: 0.5rem;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        {content}
        <div class="footer">
            <p>Generated with a Rust-based static site builder</p>
        </div>
    </div>
</body>
</html>
"#
    )
}
